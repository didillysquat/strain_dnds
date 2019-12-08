#!/usr/bin/env nextflow
// First download the .fastq files
params.sra_list = ["SRR1793320", "SRR1793321", "SRR1793322", "SRR1793323", "SRR1793324", "SRR1793325", "SRR1793326", "SRR1793327", "SRR1795737", "SRR1795735"]
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
    file "*.fasta" into ch_remove_short_iso_forms_input

    script:
    // NB that the output directory for each trinity assembly must have 'trinity' in it.
    """
    Trinity --left $corrected_read_one --right $corrected_read_two --seqType fq --max_memory 150G --CPU ${task.cpus} \\
    --min_contig_length 250 --output trinity --full_cleanup
    """
}