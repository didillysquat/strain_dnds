#!/usr/bin/env nextflow
// First download the .fastq files
params.sra_list = ["SRR1793320", "SRR1793321", "SRR1793322", "SRR1793323", "SRR1793324", "SRR1793325", "SRR1793326", "SRR1793327", "SRR1795737", "SRR1795735"]
params.bin_dir = "${workflow.launchDir}/bin"
params.launch_dir = "${workflow.launchDir}"
Channel.fromList(params.sra_list).set{ch_download_fastq}
// We are getting an error when trying to initialize all 10 of the fastq-dump requests at once.
// when running in tmux. But this doesn't seem to happen outside of tmux
process download_fastq{
    tag "${sra_id}"

    publishDir path: "raw_reads", mode: "copy"

    input:
    val sra_id from ch_download_fastq

    output:
    file "*.fastq" into ch_gzip_fastq_input

    script:
    """
    fastq-dump --defline-seq '@\$sn[_\$rn]/\$ri' --split-files ${sra_id} -O .
    """
}

// Now gzip the fastq files
process gzip_fastq{
    tag "${fastq_file_one}"

    publishDir path: "raw_reads", mode: "copy"

    input:
    tuple file(fastq_file_one), file(fastq_file_two) from ch_gzip_fastq_input

    output:
    file "*.fastq.gz" into ch_make_fastqc_pre_trim, ch_trimmomatic_input
    
    script:
    """
    gzip -f $fastq_file_one $fastq_file_two
    """
}

// Now fastqc the files
process fastqc_pre_trim{
    tag "${fastq_file}"

    publishDir path: "nf_fastqc_pre_trim", mode: "copy", saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    file fastq_file from ch_make_fastqc_pre_trim.flatten()

    output:
    file "*_fastqc.{zip,html}" into ch_fastqc_pre_trim_output

    script:
    """
    fastqc -o . $fastq_file
    """
}

// Now trim the files
process trimmomatic{
	tag "${fastq_file_one}"

	publishDir path: "nf_trimmed", mode: 'copy'

	input:
	tuple file(fastq_file_one), file(fastq_file_two) from ch_trimmomatic_input
	
	output:
	// Output that will be used for the post_trim fastqc
	// It is a flast list of all of the trimmed files
	file "*.fq.gz" into ch_fastqc_post_trim_input
	// Output that will be used for the error_correction
	// This is a list of tuples that are the 1P and 2P output files only
	tuple file("*1P.fq.gz"), file("*2P.fq.gz") into ch_rcorrect_input

	script:
	outbase = fastq_file_one.getName().replaceAll('_1.fastq.gz', '.trimmed.fq.gz')
	"""
	trimmomatic PE -threads ${params.trimmomatic_threads} -basein ${fastq_file_one} \\
		-baseout $outbase \\
		ILLUMINACLIP:${params.tru_seq_pe_fasta_path}:2:30:10:2:keepBothReads \\
		LEADING:3 TRAILING:3 MINLEN:36 HEADCROP:11
	"""
}

// Now post-trim fastqc
process fastqc_post_trim{
    tag "${fastq_file}"

    publishDir path: "nf_fastqc_post_trim", mode: "copy", saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    file fastq_file from ch_fastqc_post_trim_input.flatten()

    output:
    file "*_fastqc.{zip,html}" into ch_fastqc_post_trim_output

    script:
    """
    fastqc -o . $fastq_file
    """
}

// Now error correction
process rcorrector{
    tag "${trimmed_read_one}"
    
    cpus params.rcorrector_threads
    
    publishDir path: "nf_error_corrected", mode: "copy"
    
    input:
    tuple file(trimmed_read_one), file(trimmed_read_two) from ch_rcorrect_input

    output:
    tuple file("*1P.cor.fq.gz"), file("*2P.cor.fq.gz") into ch_trinity_input

    script:
    """
    run_rcorrector.pl -1 $trimmed_read_one -2 $trimmed_read_two -od . -t ${task.cpus}
    """
}

// Now do the trinity assembly
process trinity{
    tag "${corrected_read_one}"

    cpus params.trinity_threads
    conda "envs/nf_only_trinity.yaml"
    publishDir "nf_trinity_assembly/${corrected_read_one.getName().replaceAll('.trimmed_1P.cor.fq.gz','')}", mode: "copy"

    input:
    tuple file(corrected_read_one), file(corrected_read_two) from ch_trinity_input

    output:
    file "*.fasta" into ch_remove_short_iso_forms_trinity_input
    // Also put the corrected_read_one back into a channel and use it to derive the SRR
    // base used to make the to put the 
    val "${corrected_read_one.getName().replaceAll('.trimmed_1P.cor.fq.gz','')}" into ch_remove_short_iso_forms_srr_name_input

    script:
    // NB that the output directory for each trinity assembly must have 'trinity' in it.
    """
    Trinity --left $corrected_read_one --right $corrected_read_two --seqType fq --max_memory 150G --CPU ${task.cpus} \\
    --min_contig_length 250 --output trinity --full_cleanup
    """
}

// Now remove the short isos from the trinity assembly
// NB we were having problems getting the cache of this process to work.
// The problem was that I was using $workflow.launch_dir as the base to the /bin/remove... path
// However, when I looked at the hash logs, the value of this ($workflow.launch_dir) was giving a different hash
// hash each time. I think this was because what was getting hashed was the state of this directory.
// Obviously the state was changing very often.
// To get the hash to work in the end I have replaced the $workflow.launch_dir with a params.bin_dir
// variable that I have set at the beginning of this file using the $workflow.launch_dir variable.
// I will look to see now whether the hash logs are looking at this new variable as a string or at
// the state of this directory. If they're looking at the state of the directory then caches may
// fail when we add additional script files to the /bin directory.
process remove_short_isos{
    tag "${trinity_assembly_fasta}"

    conda "envs/nf_python_scripts.yaml"
    publishDir "nf_trinity_assembly/$srrname", mode: "copy"

    input:
    file trinity_assembly_fasta from ch_remove_short_iso_forms_trinity_input
    val srrname from ch_remove_short_iso_forms_srr_name_input

    output:
    file "*.long_iso_only.fasta" into ch_orf_prediction_input
    val srrname into ch_orf_prediction_srr_name
    
    script:
    """
    python3 ${params.bin_dir}/remove_short_isos.py $trinity_assembly_fasta
    """
}

// Do ORF prediction using transdecoder
// Concurrently rename the default names output by transdecoder
// To make the long_iso_orf files unique and related to their transcriptome we will append the srrname
process orf_prediction{
    tag "${srrname}"

    conda "envs/nf_transdecoder.yaml"
    publishDir "nf_transdecoder/$srrname", mode: "copy"

    input:
    file long_iso_trinity from ch_orf_prediction_input
    val srrname from ch_orf_prediction_srr_name

    output:
    tuple file("*.pep"), file("*.cds") into ch_remove_multi_orfs_input

    script:
    """
    TransDecoder.LongOrfs -t $long_iso_trinity -O .
    mv longest_orfs.cds ${srrname}_longest_iso_orfs.cds
    mv longest_orfs.pep ${srrname}_longest_iso_orfs.pep
    """

}

// The transdecoder output can have multiple ORFs predicted per transcript
// We will once again only keep one representative per transcript and work with this for the
// ortholog prediction
// Sonic Parnoid runs from a single directory containing all of the fastas
// To enable this we will publish each of the fastas into a single directory
process remove_multi_orfs_from_pep{
    tag "${pep_file}"
    conda "envs/nf_python_scripts.yaml"
    publishDir "nf_sonicparanoid", mode: "copy"

    input:
    tuple file(pep_file), file(cds_file) from ch_remove_multi_orfs_input

    output:
    //tuple file("*.single_orf.pep"), file(cds_file) into ch_sonic_paranoid_input
    file("*.single_orf.pep") into ch_sonicparanoid_input
    
    script:
    output_path = pep_file.getName().replaceAll("longest_iso_orfs.pep", "longest_iso_orfs.single_orf.pep")
    
    "python3 ${params.bin_dir}/unique_orfs_from_pep.py $pep_file"
}

process sonicparanoid{
    tag "sonicparanoid"
    cpus params.sonicparanoid_threads
    conda "envs/nf_sonicparanoid.yaml"
    
    input:
    // We won't actually use this input. It is just here to link
    // The processes
    file pep_file from ch_sonicparanoid_input.collect()

    output:
    file "**single-copy_groups.tsv" into ch_rename_sonicparanoid_output_input

    script:
    """
    sonicparanoid -i ${params.launch_dir}/nf_sonicparanoid -o . -t ${task.cpus}
    """
}

// The sonicparanoid file we want (single-copy_groups.tsv) is buried in a number of directories
// Here we will get rid of all of the directories and rebulish just the file
process rename_sonicparanoid_output{
    tag "rename sonicparanoid output"
    publishDir path: "nf_sonicparanoid/output_10", mode: "copy"

    input:
    file sonic_long_out from ch_rename_sonicparanoid_output_input

    output:
    file sonic_long_out into ch_screen_sonicparanoid_output_input

    script:
    "echo renamingDone"
}

// The sonic paranoid output table contains orthologs that were not found in all of the transciptomes
// We will drop these orthologs and write out the .tsv again
process screen_sonicparnoid_output{
    tag "screen sonicparanoid"
    conda "envs/nf_python_scripts.yaml"
    publishDir path: "nf_sonicparanoid/output_10", mode: "copy"

    input:
    file single_copy_tsv from ch_screen_sonicparanoid_output_input

    output:
    file "screened_orthologs.tsv" into ch_align_cds_input

    script:
    """
    python3 ${params.bin_dir}/screen_orthologs.py $single_copy_tsv
    """
}