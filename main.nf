#!/usr/bin/env/ nextflow

nextflow.enable.dsl=2

def helpMessage() {
    log.info"""
Virus Alignment Pipeline:
    
Takes an input of fastq files (zipped or unzipped), performs a de novo assembly and extracts reads that align with a given virus sequence
USAGE: nextflow run main.nf [options] --input INPUT_DIR --output OUTPUT_DIR 

OPTIONS:

--input INPUT_DIR - [Required] A directory containing paired-end fastq files

--output OUTPUT_DIR - [Required] A directory to place output files (If not existing, pipeline will create)

OPTIONAL:

    --threads INT - The number of threads that can be use to run pipeline tools in parallel

    --ref REFERENCE_FASTA - The pipeline will align contigs produced by assembly to this reference

    --minLen INT - The minimum length of a read to keep post trimming [Default = 75bp]

    --minTrimQual INT - The average basecall quality threshold below which to trim a read. During trimming, trimmomatic performs a sliding window checking the average base quality, and removing the rest of the read if it drops below this treshold. [Default = 20]
    
    """
}

// Function that checks whether a directory ends in a trailing slash.
// This is useful for directory variables that are not parsed into
// file objects in the pipeline (such as the output directory).
def checkDirectoryEnding (fileName) {
    // Grabs the last character in the directory name.
    lastChar = fileName.substring(fileName.length() - 1, fileName.length())

    // Checks whether the last character is slash
    if (lastChar != "/") {
        // If it is not a slash, add that to the directory name.
        fileName = fileName + "/"
    }
    
    // Return the directory name.
    return fileName
}

// Function creates the header for the summary file based on the parameters
// supplied by the user. Because the pipeline dynamically changes based on
// what the user specifies, the summary file must also be alter to reflect
// the analysis. This will also make incorporating new modules easier.
def createSummaryHeader (hostRef, hostIdx) {
    
    // The header will always start with the sample.
    FinalHeader = 'Sample,'

    // The summary sheet will also always contain the Raw and trimmed read counts.
    FinalHeader = FinalHeader + "Raw Reads,Trimmed Reads,Deduped Reads,"

    // Next, the user may supply host removal, which would be the next useful
    // statistic to know. Thus, if a host reference or bowtie2 index is supplied,
    // add this field to the header.
    if (hostRef != false || hostIdx != false) {
        FinalHeader = FinalHeader + "Non-Host Reads,"
    }


    // Finally, the pipeline will always report the number of contigs and scaffolds
    FinalHeader = FinalHeader + "Contigs Generated,Scaffolds Generated"
    //If a reference sequence to align to is given then the pipeline will also output # of reads that map to this target sequence
    if (params.ref != false){
        FinalHeader = FinalHeader + ",Mapped Reads"
    }

    return FinalHeader
}

// If the help parameter is supplied, link display the help message
// and quit the pipeline
params.help = false
if (params.help){
    helpMessage()
    exit 0
}


// Defines input parameters. Setting to false by default
// allows us to check that these have been set by the user.
params.input = false
params.host_fasta = false
params.host_bt2_index = false
params.ref = false
params.output = false
params.always_trim_3p_bases = 0
params.always_trim_5p_bases = 0
params.minLen = 75
params.minTrimQual = 20
params.single_read = false
params.alignmentMode = "--local"
params.phred = 33
params.mismatches_allowed = false
params.classify_singletons = false
params.max_blast_nt_evalue = "1e-10"
params.max_blastx_nr_evalue = "1e-3"
params.tally_filter_kingdom = false
//memory parameters
params.threads = 10
params.scripts_bindir = false
params.minimum_contig_length = 200




// Inports modules
include { Setup } from "./modules.nf"
include { Index_Host_Reference } from './modules.nf'
include { QC_Report } from './modules.nf'
// The same module cannot be used more than once,
// thus it is aliased to be used multiple times.
include { QC_Report as QC_Report_Trimmed } from './modules.nf'
include { QC_Report as QC_Report_Deduped } from './modules.nf'
include { QC_Report as QC_Report_Host_Removed } from './modules.nf'
include { Trimming } from './modules.nf'
include { Remove_PCR_Duplicates } from './modules.nf'
include { Host_Read_Removal } from './modules.nf'
include { Spades_Assembly } from './modules.nf'
//include { Retrieve_Contigs} from './modules.nf'
//include { Quantify_Read_Mapping} from './modules.nf'
include { Contig_Alignment } from './modules.nf'
include { Write_Summary } from './modules.nf'
adapters = file("${baseDir}/adapters.fa")

// Checks the input parameter
inDir = ""
if (params.input == false) {
    // If the parameter is not set, notify the user and exit.
    println "ERROR: No input directory provided. Pipeline requires an input directory."
    exit(1)
}
else if (!(file(params.input).isDirectory())) {
    // If the input directory is not set, notify the user and exit.
    println "ERROR: ${params.input} is not an existing directory."
    exit(1)
}
else {
    inDir = file(params.input).toString()
}

println "Input Directory: ${inDir}"

// Create a channel for the input files.

if (params.single_read != false) {
    // Single-end reads
    inputFiles_ch = Channel.fromPath("${inDir}/*.f*q*")
        .map { file -> 
            baseName = file.getBaseName() // Extract the base name of the file
            [baseName, file] // Create a tuple with the base name and the file
        }
} else {
    // Paired-end reads
    inputFiles_ch = Channel.fromFilePairs("${inDir}/*.f*q*") { file ->
        // Extract everything before the first underscore
        file.getBaseName().split("_")[0]
    }
}
// Checks the output parameter.
outDir = ''
if (params.output == false) {
    // If the parameter is not set, notify the user and exit.
    println "ERROR: No output directory provided. Pipeline requires an output directory."
    exit(1)
}
else {
    // If the parameter is set, ensure that the directory provided ends
    // in a trailing slash (to keep things consistent throughout) the
    // pipeline code.
    outDir = file(params.output).toString()
}
println outDir



// Parses the ref option.
refFile = ''
if (params.ref != false) {
    if (!(file(params.ref).exists())) {
        // If the reference file did not exist, notify the user and exit.
        println "ERROR: ${params.ref} does not exist."
        exit(1)
    }
    else {
        // Parse the provided file into a file object.
        refFile = file(params.ref)
    }
}


println "Input Directory: ${inDir}"

assembler = 'SPADES'

summaryHeader = createSummaryHeader(params.host_fasta, params.host_bt2_index)

total_deduped = 0

workflow {

    Setup( summaryHeader, params.ref, params.minLen, params.minTrimQual, outDir )

    // Use FASTQC to perform an initial QC check on the reads
    QC_Report( inputFiles_ch, outDir, "FASTQC-Pre-Processing", params.threads )

    // Perform adapter and quality trimming with Trimmomatic.
    Trimming( inputFiles_ch, outDir, adapters, params.minLen, params.minTrimQual)

    // Use FASTQC to perform a QC check on the trimmed reads.
    QC_Report_Trimmed( Trimming.out[0], outDir, "FASTQC-Trimmed", params.threads )

    // Perform PCR Duplicate removal using prinseq.
    Remove_PCR_Duplicates( Trimming.out[0], outDir, Trimming.out[1])

    // Use FASTQC to perform a QC check on the deduped reads.
    QC_Report_Deduped( Remove_PCR_Duplicates.out[0], outDir, "FASTQC-Deduplicated", params.threads )
  
    // Perform de novo assembly using spades.
    Spades_Assembly( Remove_PCR_Duplicates.out[0], outDir, params.threads, params.phred, Remove_PCR_Duplicates.out[1] )
    if (params.ref != false) {
        // Align the contigs to a reference genome using minimap2 and samtools
        Contig_Alignment( Spades_Assembly.out[0], outDir, refFile, Spades_Assembly.out[2])  
        Write_Summary( Contig_Alignment.out[1], outDir )
        }
    else{
        Write_Summary( Spades_Assembly.out[2], outDir )
    }
    
}