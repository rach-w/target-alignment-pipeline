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


    //If a reference sequence to align to is given then the pipeline will also output # of reads that map to this target sequence
    if (params.ref_fasta != false || params.ref_bt2_index != false){
        FinalHeader = FinalHeader + ",Mapped Reads,Mapped_Reads_per_Million_Unique_Reads"
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
params.ref_fasta = false
params.ref_bt2_index = false
params.output = false
params.always_trim_3p_bases = 0
params.always_trim_5p_bases = 1
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
include { Generate_Bowtie_Index as Index_Host_Reference } from './modules.nf'
include { Generate_Bowtie_Index as Index_Target_Reference} from './modules.nf'
include { QC_Report } from './modules.nf'
// The same module cannot be used more than once,
// thus it is aliased to be used multiple times.
include { QC_Report as QC_Report_Trimmed } from './modules.nf'
include { QC_Report as QC_Report_Deduped } from './modules.nf'
include { QC_Report as QC_Report_Host_Removed } from './modules.nf'
include { Trimming } from './modules.nf'
include { Remove_PCR_Duplicates } from './modules.nf'
include { Host_Read_Removal } from './modules.nf'
include { Bowtie2Alignment } from './modules.nf'
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

hostRefData = ''
hostRefIdxData = ''
hostRefName = 'NONE'
if (params.host_fasta != false && params.host_bt2_index != false) {
    // If both options are supplied, notify the user and exit.
    println "ERROR: you have specified both a host fasta file and bowtie2 index. Please only supply one."
    exit(1)
}
else if (params.host_fasta != false) {
    if (!(file(params.host_fasta).exists())) {
        // If the file supplied does not exist, notify the user and exit.
        println "ERROR: ${params.host_fasta} does not exist."
        exit(1)
    }
    else {
        // Parse the file into a file object
        hostRef = file(params.host_fasta)
        // Use the getBaseName() function to 
        // get the reference name. This will be
        // used to name the bowtie2 index.
        hostRefName = hostRef.getBaseName()
        // Place these both into a tuple.
        hostRefData = tuple(hostRefName, hostRef)
    }
}
// If the user supplied the --host_bt2_index
else if (params.host_bt2_index != false) {
    if (!(file(params.host_bt2_index).exists())) {
        // If the index provided does not exist, notify the user and exit.
        println "Error: ${params.host_bt2_index} does not exist."
        exit(1)
    }
    else {
        // Parse the directory into a file object
        hostRefDir = file(params.host_bt2_index)
        println hostRefDir
        // Grab a list of file objects from the directory
        // ending in .bt2
        indexFiles = file("${hostRefDir}/*.bt2")
        if (indexFiles.size() == 0){
            // If there are no file in the directory ending in bt2, notify the user and exit
            println "Index Directory provided (${params.host_bt2_index}) does not contain any bt2 files"
            exit(1)
        }
        else {
            // Use the getSimpleName() function to grab the base name
            // of the index files (getBaseName() removes anything following
            // the first . in a file name.)
            hostRefName = indexFiles[0].getSimpleName()
            println hostRefName
            // Place the index dir and name into a tuple.
            hostRefIdxData = tuple(hostRefDir, hostRefName)
        }
    }
}

// Parses the ref option.
targetRefData = ''
targetRefIdxData = ''
targetRefName = 'NONE'
if (params.ref_fasta != false && params.ref_bt2_index != false) {
    // If both options are supplied, notify the user and exit.
    println "ERROR: you have specified both a host fasta file and bowtie2 index. Please only supply one."
    exit(1)
}
else if (params.ref_fasta != false) {
    if (!(file(params.ref_fasta).exists())) {
        // If the file supplied does not exist, notify the user and exit.
        println "ERROR: ${params.ref_fasta} does not exist."
        exit(1)
    }
    else {
        // Parse the file into a file object
        targetRef = file(params.ref_fasta)
        // Use the getBaseName() function to 
        // get the reference name. This will be
        // used to name the bowtie2 index.
        targetRefName = targetRef.getBaseName()
        // Place these both into a tuple.
        targetRefData = tuple(targetRefName, targetRef)
    }
}
// If the user supplied the --host_bt2_index
else if (params.target_bt2_index != false) {
    if (!(file(params.target_bt2_index).exists())) {
        // If the index provided does not exist, notify the user and exit.
        println "Error: ${params.target_bt2_index} does not exist."
        exit(1)
    }
    else {
        // Parse the directory into a file object
        targetRefDir = file(params.target_bt2_index)
        println targetRefDir
        // Grab a list of file objects from the directory
        // ending in .bt2
        indexFiles = file("${targetRefDir}/*.bt2")
        if (indexFiles.size() == 0){
            // If there are no file in the directory ending in bt2, notify the user and exit
            println "Index Directory provided (${params.target_bt2_index}) does not contain any bt2 files"
            exit(1)
        }
        else {
            // Use the getSimpleName() function to grab the base name
            // of the index files (getBaseName() removes anything following
            // the first . in a file name.)
            targetRefName = indexFiles[0].getSimpleName()
            println targetRefName
            // Place the index dir and name into a tuple.
            targetRefIdxData = tuple(targetRefDir, targetRefName)
        }
    }
}


println "Input Directory: ${inDir}"

assembler = 'SPADES'

summaryHeader = createSummaryHeader(params.host_fasta, params.host_bt2_index)

total_deduped = 0

workflow {

    Setup( summaryHeader, params.ref_fasta, params.minLen, params.minTrimQual, outDir )

    if (params.host_fasta){
        Index_Host_Reference(hostRefData, outDir, params.threads )
        Index_Host_Reference.out[0]
            .set{hostRefIdxData}
    }

    if (params.ref_fasta){
        Index_Target_Reference(targetRefData, outDir, params.threads)
        Index_Target_Reference.out[0]
            .set{targetRefIdxData}
    }
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

    //Remove host reads
    Host_Read_Removal(Remove_PCR_Duplicates.out[0], outDir, hostRefIdxData, params.alignmentMode, params.threads, Remove_PCR_Duplicates.out[1])

    // Align the contigs to a reference bowtie2 index
    Bowtie2Alignment( Host_Read_Removal.out[0], outDir, targetRefIdxData, params.alignmentMode, params.threads, Host_Read_Removal.out[2], Remove_PCR_Duplicates.out[2])  

    //Write to the summary file with stats
    Write_Summary( Bowtie2Alignment.out[1], outDir )
        
    
}