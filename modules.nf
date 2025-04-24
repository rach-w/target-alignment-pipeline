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
process Index_Host_Reference {
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
    //The minimum read lenths to be allowed post trimming
    val minLen
    // The threshold below which a read is trimmed while
    // performing sliding window trimming
    val minTrimQual 

    output:
    tuple val(base), path("*_f.fastq") 
    tuple val(base), path("*_f.fastq")
    env 'summary' 

    publishDir "${outDir}/trimming", mode: "copy"
    // TODO: parameterize adapter sequences
    script:

    // this handles paired-end data, in which case must specify a paired output file
    def paired_output   = initial_fastq[1] ? "-p ${base}_R2_f.fastq" : ""
    def paired_adapters = initial_fastq[1] ? "-A AGATCGGAAGAGC -G GCTCTTCCGATCT -A AGATGTGTATAAGAGACAG -G CTGTCTCTTATACACATCT" : ""
    // TODO: don't trim this much for non-amplicon data!
    def paired_trimming = initial_fastq[1] ? "-U $params.always_trim_5p_bases -U -${params.always_trim_3p_bases}" : ""

    """
    #!/bin/bash
    raw_reads_1=\$((\$(zcat -f ${initial_fastq} | wc -l)/4))
    raw_reads_2=\$((\$(zcat -f ${initial_fastq} | wc -l)/4))

    total_raw=\$((\$raw_reads_1 + \$raw_reads_2))
    
    cutadapt \
    -a AGATCGGAAGAGC -g GCTCTTCCGATCT -a AGATGTGTATAAGAGACAG -g CTGTCTCTTATACACATCT \
    $paired_adapters \
    -q 30,30 \
    --minimum-length ${params.post_trim_min_length} \
    -u ${params.always_trim_5p_bases} \
    -u -${params.always_trim_3p_bases} \
    $paired_trimming \
    -o ${base}_R1_f.fastq \
    $paired_output \
    $initial_fastq 
    --info-file trim_summary.tsv
    
    trimmed_reads_1=\$((\$(wc -l < ${base}_R1_f.fastq)/4))
    trimmed_reads_2=\$((\$(wc -l < ${base}_R2_f.fastq)/4))

    total_trimmed="\$((\$trimmed_reads_1 + \$trimmed_reads_2))"

    summary="${base}, \$total_raw, \$total_trimmed"
    

    """

}
// Removes PCR Duplicates from a set of reads using prinseq.
process Remove_PCR_Duplicates {
    input:
    //tuple with sample id/basename & fastq
    tuple val(base), path(input_fastq)
    // The output directory
    val outDir
    // The existing summary file
    val existingSummary

    val total_deduped

    output:
    //tuple with base & dedupelicated fastq file
    tuple val(base), path("*_fu.fastq") 
    //summary string to be written
    env 'summary'
    //number of unique reads to be used for reads mapped per million unique reads calculation
    env 'total_deduped'

    script:

    // this handles paired-end data, in which case must specify a paired output file
    def r1 = input_fastq[0]
    def r2 = input_fastq[1] 
    def prefix_param    = input_fastq[1] ? "-u 30" : "-u 30"
    def paired_input    = input_fastq[1] ? "-i2 $r2" : ""
    def paired_output   = input_fastq[1] ? "-o2 ${base}_R2_fu.fastq" : ""

    """
    #!/bin/bash

    cd-hit-dup \
    -e $params.mismatches_allowed \
    $prefix_param \
    -i $r1 \
    $paired_input \
    -o ${base}_R1_fu.fastq \
    $paired_output \

    deduped_reads_1=\$((\$(wc -l < ${base}_R1_fu.fastq ) /4 ))
    deduped_reads_2=\$((\$(wc -l < ${base}_R2_fu.fastq) /4 ))

    total_deduped=\$(( \$deduped_reads_1 + \$deduped_reads_2 ))
    summary="${existingSummary}, \$total_deduped"

    """
}

// Aligns reads to a reference gneome using bowtie2
process Bowtie2Alignment {
    input:
        // Tuple contains the file basename and paired-end reads
        tuple val(base), path(trimmed_fastq)
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

    output:
        // Tuple contains the file basename and the alignment in a sorted bam file
        tuple val(base), file("${base}.bam")
        // The summary string containing the number of mapped reads
        env 'summary'
    
    publishDir "${outDir}/${base}-Intermediate-Files/", mode: 'copy', pattern: "${base}.bam"

    script:
    /*
    Aligns the reads to the reference genome using bowtie2. ocal alignment is used (--local)
    to ensure the reads are aligne without gaps.

    The alignment is then converted to bam format and sorted using samtools.

    Samtools is then used to grab the number of mapped reads from the alignment,
    and this is added to the summary string.
    */
    """
    #!/bin/bash

    bowtie2 --threads ${threads} -x ${refDir}/${refName} -1 ${R1} -2 ${R2} ${alignmentMode} -S ${base}-align.sam

    samtools view -b ${base}-align.sam | samtools sort > ${base}.bam

    mapped_reads="\$(samtools view -F 0x04 -c ${base}.bam)"

    summary="${existingSummary}, \$mapped_reads"
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
    def r1 = input_fastq[0] 
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

    mv ${base}_host_removed.1 ${base}_host_removed_1.fastq
    mv ${base}_host_removed.2 ${base}_host_removed_2.fastq

    nonHost_reads_1=\$((\$(wc -l < ${base}_host_removed_1.fastq)/4))
    nonHost_reads_2=\$((\$(wc -l < ${base}_host_removed_2.fastq)/4))

    total_nonHost=\$((\$nonHost_reads_1 + \$nonHost_reads_2))

    summary="${existingSummary}, \$total_nonHost"
    
    """
}
// Uses spades to produce a de novo assembly.
process Spades_Assembly {
    input:
        // Tuple contains the file basename and paired-end reads.
        tuple val(base), path(input_fastq)
        // The output directory
        val outDir
        // The number of threads provided.
        val threads
        //phred offset parameter for spades
        val phred
        // The existing summary string
        val existingSummary
    output:
        // Tuple contains the basename of the sample and the assembled contigs 
        // produced by spades.
        tuple val(base), file("${base}-contigs.fasta")
        // The scaffolds produced by spades
        file "${base}-scaffolds.fasta"
        // The summary string containing the number of contigs and scaffolds
        env 'summary'

        tuple val(base), path(input_fastq)
    
    publishDir "${outDir}/assembled_contigs", mode: 'copy'

    script:
    /*
    Runs spades using the provided paired-end reads.

    The contigs and scaffolds are renamed and moved, and
    the number of contigs/scaffolds are recorded.
    
    The number of contigs and scaffolds are added to the summary string.
    */
    // handle single-end or paired-end inputs
    def r1 = input_fastq[0] 
    def r2 = input_fastq[1] ? input_fastq[1] : ""
    def spades_input  =      input_fastq[1] ? "-1 $r1 -2 $r2" : "-s $r1"
    """
    #!/bin/bash

    spades.py --threads ${threads} ${spades_input} -o ${base}-Assembly --phred-offset ${phred}

    if [[ -f "${base}-Assembly/contigs.fasta" ]]; then
        mv ${base}-Assembly/contigs.fasta ./${base}-contigs.fasta
    else
        touch ./${base}-contigs.fasta
    fi

    num_contigs=\$(grep ">" ${base}-contigs.fasta | wc -l)

    if [[ -f "${base}-Assembly/scaffolds.fasta" ]]; then
        mv ${base}-Assembly/scaffolds.fasta ./${base}-scaffolds.fasta
    else
        touch ./${base}-scaffolds.fasta
    fi

    num_scaffolds=\$(grep ">" ${base}-scaffolds.fasta | wc -l)

    summary="${existingSummary},\$num_contigs, \$num_scaffolds"
    """
}


// Get contigs from assembly, converting into 1-line fasta format
/*
process Retrieve_Contigs {
  
  label 'lowmem'  

  input:
  tuple val(base), file("${base}-contigs.fasta")
  tuple val(base), path(input_fastq)
  val outDir 

  output:
  tuple val(base), path("${base}_contigs.fa")
  tuple val(base), path(input_fastq) 

  publishDir "${outDir}/contig_oneline", mode:'link', pattern:"*.fa"

  script:
  """
  #!/bin/bash

  # consolidate assembly output
  # this forces contigs fasta into 1-line format
  # also remove contigs shorter than minimum_contig_length bases
  seqtk seq -A -L ${params.minimum_contig_length} ${base}-contigs.fasta > ${base}_contigs.fa
  """
}
*/

/*
  Map reads to contigs so that contigs can be weighted by the # of 
  reads they account for.
*/
/*
process Quantify_Read_Mapping {
  label 'lowmem_threaded'
  

  input:
  //basename & path to one-line contigs output by read mapping
  tuple val(base), path(contigs)
  //basename & path to input fastq file
  tuple val(base), path(input_fastq) 
  //number of threads to use 
  val (threads)
  //directory to output to
  val outDir

  output:
  tuple val(base), path("${base}_contigs.fasta")
  path("${base}_contig_weights.txt")
  path("contig_mapping.sam")

  publishDir "${outDir}/read_mapping", mode:'copy'

  script:
  def r1 = input_fastq[0] 
   
  """
  #!/bin/bash

  bowtie2-build ${contigs} contig_bt_index

  # C,120,1 makes min score 120 -> ~corresponds to 100% identity over ~60 bases
  bowtie2 -x contig_bt_index --local --score-min C,120,1 -q -U ${r1} -p ${threads} -S contig_mapping.sam 

  # count the # of reads that hit each contig 
  

  awk '( \$0 !~ /^@/ ) && ( \$3 != "*" ) { print }' contig_mapping.sam |   # ignore header lines (starting w/ @) and unmapped reads (field 3 == *)
  cut -f 3 |                 # cut out 3rd field (subject ID)
  sort |                     # sort subject IDs
  uniq -c |                  # collapse to unique list of subject IDs and tally counts
  sort -nr |                 # sort numerically by tally, largest first
  awk '{print \$2 "\t" \$1}' > ${base}_contig_weights.txt   # reverse column order and write to output file

  mv ${base}_contigs.fa  ${base}_contigs.fasta
  """
}
*/


// Generates an alignment of the asembly contigs to a reference genome *if provided
// using minimap.
process Contig_Alignment {
    input:
        // Tuple contains the sample basename
        // and the assembly fasta to align
        tuple val(base), file(assembly)
        // The output directory
        val outDir
        // the reference fasta file to be aligned to.
        file ref

    output:
        // Tuple contains the file basename and the alignment bam file.
        tuple val(base), file("${base}-contig-align.bam")

    publishDir "${outDir}", mode: 'copy'
    
    script:
    /*
    Uses minimap2 to align the contigs to the reference fasta.

    Then uses samtools to conver the alignment sam into bam format
    (samtools view). The bam format is then sorted and stored in a bam
    file (samtools sort).
    */
    """
    #!/bin/bash

    minimap2 -ax asm5 ${ref} ${assembly} > align.sam

    samtools view -b align.sam | samtools sort > ${base}-contig-align.bam

    
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




/*
   Output a matrix of # of virus-mapping reads
*/
/*
process Virus_Mapping_Matrix {
  label 'lowmem_nonthreaded'

  input:
  path(normalized_tally_files) 
  val(outDir)

  output:
  path("*.tsv") 
  publishDir "${outDir}/virus_mapping_matrix", mode:'copy'

  script:
  """
  Rscript ${params.scripts_bindir}/taxa_matrix.R ${normalized_tally_files} -v
  """
}

process Create_Heatmap{
  label 'lowmem_nonthreaded'

  input:
  path(mapping_matrix)
  val(outDir)

  output:
  path("*.pdf")
  publishDir "${outDir}/heatmap", mode: 'copy'

  script:
  """
  Rscript ${params.scripts_bindir}/generate_heatmap.R -i ${mapping_matrix} 
  """
}
*/






