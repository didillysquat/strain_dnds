#!/usr/bin/env nextflow
params.bin_dir = "${workflow.launchDir}/bin"
params.launch_dir = "${workflow.launchDir}"
Channel.fromList(params.sra_list).set{ch_download_fastq}
// We are getting an error when trying to initialize all 10 of the fastq-dump requests at once.
// when running in tmux. But this doesn't seem to happen outside of tmux
process download_fastq{
    tag "${sra_id}"
    conda "envs/nf_general.yaml"
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
    conda "envs/nf_general.yaml"
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
    conda "envs/nf_general.yaml"
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
    conda "envs/nf_general.yaml"
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
    conda "envs/nf_general.yaml"
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
    conda "envs/nf_general.yaml"
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
    conda "envs/nf_general.yaml"
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
    file cds_file into ch_write_unaligned_cds_fastas_fas_input
    
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
    file "screened_orthologs.tsv" into ch_write_unaligned_cds_fastas_tab_input

    script:
    """
    python3 ${params.bin_dir}/screen_orthologs.py $single_copy_tsv
    """
}


// Work through the screened_orthologs and for each ortholog make a directory
// and write out a fasta that contains the sequence for each of the strains
// Input to the script will be the .tsv. Each of the .cds fastas will be
// in the nextflow working directory due to the channel input
// The column titles of the screened ortholog .tsv are the .pep file
// names (e.g. SRR1793320_longest_iso_orfs.single_orf.pep)
// The cds files that we will be getting the sequences from are all prefixed
// with the SRRXXX (e.g. SRR1793320_longest_iso_orfs.cds). As such we can use
// this SRRXXX in the python script to map between the table and the .cds files
process write_unaligned_cds_fastas{
    tag "write_unaligned_cds_fastas"
    conda "envs/nf_python_scripts.yaml"
    publishDir path: "nf_local_alignments", mode: "copy"
    
    input:
    file screened_orth_table from ch_write_unaligned_cds_fastas_tab_input
    file cds_fastas from ch_write_unaligned_cds_fastas_fas_input.collect()

    output:
    file "**/*_unaligned_cds.fasta" into ch_align_using_guidance_input

    script:
    """
    python3 ${params.bin_dir}/write_out_unaligned_cds_fastas.py $screened_orth_table
    """
}

// Align each of the unaligned cds files.
// We will do this in two processes.
// First we will perform the guidance analysis
// Then we will use the three output files to create the aligned fastas
// That have been cropped and had the appropriate low quality columns dropped
// As always, Guidance is proving tricky to get to run
// We have had to explicityly add per-bioperl to the env yaml
// Then, the bash wrapper for the perl execution of the guidance.pl
// is not working so you have to explicityly run perl and the full path to guidance.pl
// Because we can't know what the full path is, we will write a quick python script
// to find this out and run it.
// Guidance seems to be quite unstable. It seems to produce errors for no reason.
// However, the pipeline can simply be restarted and everything seems to behave OK.
// As such I have now implemented an automatic retry process where the python
// script will attempt to re-run the guidance analysis in question before quitting.
process align_using_guidance{
    tag "${unaligned_fasta.toString().split('_')[0]}"
    conda "envs/nf_guidance.yaml"

    input:
    file unaligned_fasta from ch_align_using_guidance_input.flatten()

    output:
    tuple file("*.MAFFT.Guidance2_res_pair_res.PROT.scr"), file("*.MAFFT.PROT.aln"), file("*.MAFFT.aln.With_Names") into ch_process_guidance_output_input

    script:
    orth_group_id = unaligned_fasta.toString().split('_')[0]
    """
    python3 ${params.bin_dir}/run_guidance.py $unaligned_fasta
    """
}

// This process will use the three ouput files from align_using_guidance
process process_guidance_output{
    // We want to publish each of the output files into the corresponding ortholog directory
    // So here we will attempt to use the closure to extract the ortholog number
    /// and assign the publication directory accrodingly.
    tag "${aa_cols_score_file_path.toString().split('_')[0]}"
    publishDir path: "nf_local_alignments", mode: "copy", saveAs:   {filename -> def orth_id = filename.split('_')[0]
                                                                                return "$orth_id/$filename"}

    input:
    tuple file(aa_cols_score_file_path), file(aa_alignment_file_path), file(cds_alignment_file_path) from ch_process_guidance_output_input

    output:
    file "*_cropped_aligned_aa.fasta" into ch_model_test_input
    file "*_cropped_aligned_cds.fasta" into ch_run_codeml_align_input

    script:
    """
    python3 ${params.bin_dir}/process_guidance_output.py $aa_cols_score_file_path $aa_alignment_file_path $cds_alignment_file_path
    """
}

// Now find the best protein evo model
// In some cases the fasta files that resulted from the process_guidance_output may be empty
// As such we will wrap the modeltest-ng running in a python script that will check for this.
// If there is this problem with the aligned fasta then we will return code 0
// before making sure to delete any *.out file that might have been made
// We will make it so that the *_prottest_result.out is output empty.
// We can then check for the empty .out file in the make_master_alignment
process model_test{
    tag "${cropped_aligned_aa_fasta.toString().split('_')[0]}"
    conda "envs/nf_modeltest-ng.yaml"
    publishDir path: "nf_prot_out"

    input:
    file cropped_aligned_aa_fasta from ch_model_test_input

    output:
    tuple file("*_prottest_result.out"), file(cropped_aligned_aa_fasta) into ch_make_master_alignment_input
    val "nf_prot_out" into ch_make_master_alignment_dir_input
    script:
    """
    python3 ${params.bin_dir}/run_model_test.py $cropped_aligned_aa_fasta
    """
}

// ch_make_master_alignment_input.collect().subscribe {  println "Got: $it"  }

// To make the master tree we will work witha single process that
// will need to iterthrough each of the protein model outputs.
// It will also need access to the aa cropped and alignment files
// We will supply both of these in two seperate input channels
// The output will be a master fasta and a q file for raxml that delimits the partions
// that can then be fed into the treemaking
process make_master_alignment_and_q_file{
    tag "make_master_alignment"
    conda "envs/nf_python_scripts.yaml"

    input:
    val nf_prot_out_dir from ch_make_master_alignment_dir_input.collect()

    output:
    tuple file("master_fasta_for_tree.fasta"), file("q_partition_file.q") into ch_make_tree_input

    script:
    """
    python3 ${params.bin_dir}/make_master_alignment.py ${params.launch_dir}/${nf_prot_out_dir[0]}
    """
}

// Make a ML tree using raxml
process make_tree{
    tag "make_tree"
    cpus params.raxml_threads
    conda "envs/nf_raxml.yaml"

    input:
    tuple file(master_fasta), file(q_partition_file) from ch_make_tree_input

    output:
    tuple file("RAxML_bestTree.strain_dn_ds"), file("RAxML_bipartitionsBranchLabels.strain_dn_ds"), file("RAxML_bipartitions.strain_dn_ds") into ch_annotate_tree_input
    file "RAxML_bestTree.strain_dn_ds" into ch_run_codeml_tree_input
    script:
    """
    raxmlHPC-PTHREADS-AVX2 -s $master_fasta -q $q_partition_file -x 183746 -f a, -p \\
    83746273 -# 1000 -T ${task.cpus} -n strain_dn_ds -m PROTGAMMAWAG
    """
}

// Apply more human readable labels to the tree
// Currently the labels are the SRRXXX this will
//add the strain and species to the labels
process annotate_tree{
    tag "annotate_tree"
    conda "envs/nf_python_scripts.yaml"
    publishDir path: "nf_master_tree"

    input:
    tuple file(tree_one), file(tree_two), file(tree_three) from ch_annotate_tree_input

    output:
    tuple file("RAxML_bestTree.strain_dn_ds_named"), file("RAxML_bipartitionsBranchLabels.strain_dn_ds_named"), file("RAxML_bipartitions.strain_dn_ds_named") into ch_tree_for_dnds_input

    script:
    """
    python3 ${params.bin_dir}/annotate_tree.py $tree_one $tree_two $tree_three
    """
}

// ch_run_codeml_align_input.combine(ch_run_codeml_tree_input).subscribe {  println "Got: $it"  }
// Here we create codeml control files.
// Previously we were going to a lot of effort to incorporate
// some sort of parallelisation. However, let's see if we can make use of nextflow
// format here and run one instance per cds input.
// When we run this we will check to see that there are sequences in the alignment and 
// that the alignment is divisible by 3. If either of these assertions fails
// we will exit without error. This means that the *.out file output needs to be optional
process run_codeml{
    tag "$cds_file"
    conda "envs/nf_codeml.yaml"
    publishDir path: "nf_codeml_out"
    input:
    tuple file(cds_file), file(tree) from ch_run_codeml_align_input.combine(ch_run_codeml_tree_input)
    
    output:
    file "*_codeml_results.out" optional true into ch_collate_codeml_results_intput
    script:
    """
    python3 ${params.bin_dir}/run_codeml.py $cds_file $tree
    """
}

// ch_collate_codeml_results_intput.collect().subscribe {println "Got: $it"}
process collate_codeml_results{
    tag "collate_codeml_results"
    conda "envs/nf_python_scripts.yaml"

    input:
    file codeml_out_file from ch_collate_codeml_results_intput.collect()
    output:
    file "codeml_results_df.csv" into ch_collate_codeml_results

    script:
    """
    python3 ${params.bin_dir}/collate_codeml_results.py ${params.sra_list_as_csv}
    """
}