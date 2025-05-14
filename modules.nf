// Creates a parameters file and a summary file to 
// be added to later
process Setup {
    input:
        // The header to write to the summary file.
        val summaryHeader
        // The name of the reference supplied (for use
        // in the parameters file)
        val refName
        // The minimum read length allowed post trimmming (for
        // use in the parameters file)
        val minLen
        // The threshold below which a read is trimmed while
        // performing sliding window trimming
        val minTrimQual
        // The output directory to be used.
        val outDir
        
    output:
        // The parameters file created.
        file "analysis-parameters.txt"
        // The blank summary file to be added to.
        file "stats-summary.csv"

    publishDir "${outDir}", mode: 'copy'

    script:
    /*
    Creates a parameters file (in case the user needs to look back at how they ran the analysis)
    as well as a blank summary file, which will later be populated with data from each
    sample.

    The parameters file contains:
        1. The name of the reference supplied
        2. The minimum read length allowed after trimmming
        3. The threshold below which a read is trimmed 

    The summary file contains:
        1. The sample
        2. Raw Reads
        3. Reads after trimming
        4. Reads after deduplication
        5. Reads after host removal
        6. Number of contigs
        7. Number of scaffolds
        8. Number of reads that mapped to the target sequence
    */
    """
    #!/bin/bash

    touch analysis-parameters.txt

    echo "Minimum Read Length Allowed : ${minLen} bp" >> analysis-parameters.txt
    echo "Trimming Quality Threshold : ${minTrimQual}" >> analysis-parameters.txt
    echo "Target genome : ${refName}" >> analysis-parameters.txt

    touch stats-summary.csv

    echo "${summaryHeader}" > stats-summary.csv
    """
}

// Builds a bowtie2 index for a provided reference file
process Generate_Bowtie_Index {
    input:
        // Tuple contains the reference name and reference file.
        tuple val(refName), file(ref)
        // The output directory
        val outDir
        // The number of threads provided
        val threads

    output:
        // Tuple containing the reference index directory and the index file name.
        tuple file("ref-idx/"), val(refName)

    publishDir "${outDir}", mode: 'copy'

    /*
    Creates an index Directory and then runs bowtie2-build to create the index
    */
    script:
    """
    #!/bin/bash
    
    mkdir ref-idx/

    bowtie2-build --threads ${threads} -f ${ref} ref-idx/${refName}
    """
}

// Creates a fastqc report for a set of reads provided.
process QC_Report {
    input:
        // Tuple contains the file basename as well as the read files
        tuple val(base), path (fastq)
        // The output directory
        val outDir
        // The name of the directory to place the fastqc output into (allows
        // for this command to be used multiple times and to separate the output
        // i.e. pre-processed-reads vs trimmed-reads.)
        val dirName
        // The number of threads provided.
        val threads

    output:
        // The directory cotaining the fastqc files produced.
        path("${base}")
    
    publishDir "${outDir}/${dirName}/", mode: 'copy'

    // Creates a Directory and then runs fastqc on a set of reads.
    script:
    """
    #!/bin/bash

    mkdir ${base}

    fastqc $fastq -o ${base} --threads ${threads}
    """
}


// Performs quality and adapter trimming on a set of paired-end reads.
process Trimming {
    input:
    tuple val(base), path(initial_fastq)
    // The output directory
    val outDir
    // The adapter file in fasta format.
    file adapters
    //The minimum read lengths to be allowed post trimming
    val minLen
    // The threshold below which a read is trimmed while
    // performing sliding window trimming
    val minTrimQual 

    output:
    tuple val(base), path("*_f.fastq"), emit: trimmed_reads
    env 'summary', emit: summary

    publishDir "${outDir}/trimming", mode: "copy"

    script:
    def is_paired = initial_fastq instanceof List && initial_fastq.size() > 1
    
    // Handle paired-end specific parameters
    def paired_output = is_paired ? "-p ${base}_R2_f.fastq" : ""
    def paired_adapters = is_paired ? "-A AGATCGGAAGAGC -G GCTCTTCCGATCT -A AGATGTGTATAAGAGACAG -G CTGTCTCTTATACACATCT" : ""
    def paired_trimming = is_paired ? "-U $params.always_trim_5p_bases -U -${params.always_trim_3p_bases}" : ""
    
    def input_file = is_paired ? initial_fastq[0] : initial_fastq

    """
    #!/bin/bash
    raw_reads=\$((\$(zcat -f ${input_file} | wc -l)/4))
    
    cutadapt \
    -a AGATCGGAAGAGC -g GCTCTTCCGATCT -a AGATGTGTATAAGAGACAG -g CTGTCTCTTATACACATCT \
    $paired_adapters \
    -q 30,30 \
    --minimum-length ${params.minLen} \
    -u ${params.always_trim_5p_bases} \
    -u -${params.always_trim_3p_bases} \
    $paired_trimming \
    -o ${base}_R1_f.fastq \
    $paired_output \
    $input_file \
    --info-file trim_summary.tsv
    
    trimmed_reads=\$((\$(wc -l < ${base}_R1_f.fastq)/4))
    
    summary="${base}, \$raw_reads, \$trimmed_reads"
    """
}

process Remove_PCR_Duplicates {
    input:
    // tuple with sample id/basename & fastq
    tuple val(base), path(input_fastq)
    // The output directory
    val outDir
    // The existing summary file
    val existingSummary

    output:
    // tuple with base & deduplicated fastq file
    tuple val(base), path("*_fu.fastq"), emit: deduped_reads
    // summary string to be written
    env 'summary', emit: summary
    // number of unique reads to be used for reads mapped per million unique reads calculation
    env 'total_deduped', emit: total_deduped

    script:
    def is_paired = input_fastq instanceof List && input_fastq.size() > 1
    
    def r1 = is_paired ? input_fastq[0] : input_fastq
    def prefix_param = "-u 30"  // Same for both single and paired
    def paired_input = is_paired ? "-i2 ${input_fastq[1]}" : ""
    def paired_output = is_paired ? "-o2 ${base}_R2_fu.fastq" : ""

    """
    #!/bin/bash

    cd-hit-dup \
    -e ${params.mismatches_allowed} \
    $prefix_param \
    -i $r1 \
    $paired_input \
    -o ${base}_R1_fu.fastq \
    $paired_output \

    deduped_reads=\$((\$(wc -l < ${base}_R1_fu.fastq ) /4 ))
    
    if  $is_paired
    then    
        deduped_reads_2=\$((\$(wc -l < ${base}_R2_fu.fastq) /4 ))
        total_deduped=\$(( \$deduped_reads + \$deduped_reads_2 ))
        summary="${existingSummary}, \$total_deduped"
    else    
        total_deduped=\$deduped_reads
        summary="${existingSummary}, \$total_deduped"
    fi
    """
    
}


// Removes host reads by aligning to a host genome using bowtie2.
process Host_Read_Removal {
    input:
        // Tuple contains the file basename and paired-end reads
        tuple val(base), path(input_fastq) 
        // The output directory name
        val outDir
        // Tuple contains the bt2 index directory and basename of the index files.
        tuple path(refDir), val(refName) 
        //the alignment mode parameter
        val alignmentMode
        // The number of threads provided.
        val threads
        // The existing statistics string to be added to.
        val existingSummary
        
    output:
        // Tuple contains the file basename and the paired-end read files with host reads removed.
        tuple val(base), path("${base}_*.fastq")
        // A directory containing the alignment to the host file.
        file "host-reads/${base}-host.sam"
        // The summary string containing the number of reads post host removal
        env 'summary'
    
    publishDir "${outDir}/${base}-Intermediate-Files/Processed-Reads", mode: 'copy'

    script:
    /*
    Aligns the reads to a host reference genome using bowtie2. Local alignment is used (--local)
    to ensure the reads are aligne without gaps. The unaligned reads are sent back ot into a
    fastq files (--un-conc).

    The unaligned read files are then renamed to give them the typical paired end
    read name scheme.

    Finally, the number of forward and reverse reads are grabbed and used to calculate
    the total number of reads post host removal. This value is added to the summary string.
    */

    // handle single-end or paired-end inputs
    def is_paired = input_fastq instanceof List && input_fastq.size() > 1
    def r1 = is_paired ? input_fastq[0] : input_fastq
    def r2 = input_fastq[1] ? input_fastq[1] : ""
    def bowtie_file_input  = input_fastq[1] ? "-1 $r1 -2 $r2" : "-U $r1"
    def bowtie_file_output = input_fastq[1] ? "--un-conc ${base}_host_removed" : "--un ${base}_host_removed"

    """
    #!/bin/bash
    mkdir host-reads
    
    bowtie2 -x ${refDir}/${refName} \
    -q $bowtie_file_input \
    --threads ${threads}\
    $bowtie_file_output \
    --local -S host-reads/${base}-host.sam

    if  $is_paired
    then    
        mv ${base}_host_removed.1 ${base}_host_removed_1.fastq
        mv ${base}_host_removed.2 ${base}_host_removed_2.fastq

        nonHost_reads_1=\$((\$(wc -l < ${base}_host_removed_1.fastq)/4))
        nonHost_reads_2=\$((\$(wc -l < ${base}_host_removed_2.fastq)/4))

        total_nonHost=\$((\$nonHost_reads_1 + \$nonHost_reads_2))

        summary="${existingSummary}, \$total_nonHost"
    else    
        mv ${base}_host_removed ${base}_host_removed_1.fastq
        nonHost_reads=\$((\$(wc -l < ${base}_host_removed_1.fastq)/4))
        total_nonHost=\$nonHost_reads
        summary="${existingSummary}, \$total_nonHost"
    fi
    """
}

// Aligns reads to a reference gneome using bowtie2
process Bowtie2Alignment {
    label 'lowmem_threaded'
    input:
        // Tuple contains the file basename and paired-end reads
        tuple val(base), path(input_fastq)
        // The output directory name
        val outDir
        // Tuple contains the bowtie2 index directory and the name of the reference used
        tuple file(refDir), val(refName)
        // The alignment mode parameter
        val alignmentMode
        // The number of threads provided.
        val threads
        // The existing statistics string to be added to.
        val existingSummary
        // Total number of deduped reads for reads/million unique reads calculation
        val total_deduped

    output:
        // Tuple contains the file basename and the alignment in a sorted bam file
        tuple val(base), file("${base}.bam")
        // The summary string containing the number of mapped reads
        env 'summary'
    
    publishDir "${outDir}/${base}-Intermediate-Files/", mode: 'copy', pattern: "${base}.bam"
    
    script:
    //Adjusting input and output to handle paired end and single end data
    def r1 = input_fastq[0] 
    def r2 = input_fastq[1] ? input_fastq[1] : ""
    def bowtie_file_input  = input_fastq[1] ? "-1 $r1 -2 $r2" : "-U $r1"

    /*
    Aligns the reads to the reference genome using bowtie2. ocal alignment is used (--local)
    to ensure the reads are aligne without gaps.

    The alignment is then converted to bam format and sorted using samtools.

    Samtools is then used to grab the number of mapped reads from the alignment,
    and this is added to the summary string.
    */
    """
    #!/bin/bash

    bowtie2 --threads ${threads} -x ${refDir}/${refName} -q $bowtie_file_input ${alignmentMode} -S ${base}-align.sam

    samtools view -b ${base}-align.sam | samtools sort > ${base}.bam

    mapped_reads="\$(samtools view -F 0x04 -c ${base}.bam)"

    reads_per_mill=\$(echo "scale=4; (\$mapped_reads / ${total_deduped}) * 1000000" | bc)

    summary="${existingSummary}, \$mapped_reads, \$reads_per_mill"
    """
}


// Writes a line to the summary file for the sample.
process Write_Summary {
    input:
        // The summary string containing statistics collected as
        // the pipeline ran.
        val summary
        // The output directory.
        val outDir

    script:
    /*
    The summary statistics are written to the summary file.
    */
    """
    #!/bin/bash

    echo "${summary}" >> ${outDir}/stats-summary.csv
    """  

}










