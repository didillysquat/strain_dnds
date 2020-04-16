#!/usr/bin/env nextflow
//TODO for first test we will try to get 4 + 2 working.

// I have created directories that hold symlinks to either the 4 seq set or to the 10 seq set.
// Use these directories to work with either the 4 or 10 samples set.
// For the 4 or 4 + 2 analysis:
params.raw_reads_dir = "/home/humebc/projects/parky/breviolum_transcriptomes/raw_reads/4"
// For the 10 or 10 + 2 analysis:
// params.raw_reads_dir = "/home/humebc/projects/parky/breviolum_transcriptomes/raw_reads/10"
// params.sra_list = ["b_min_c", "b_min_c", "b_aenig", "b_pseudo"]
// For the 10 or 10 + 2 analysis:
// params.sra_list = ["b_min_0", "b_min_1", "b_min_2", "b_min_3", "b_min_0", "b_min_1", "b_min_2", "b_min_3", "b_aenig", "b_pseudo"]
// This is a bool of whether we are adding the b_5 and the b_fav transcriptomes as well
params.additional_transcriptomes = true
params.bin_dir = "${workflow.launchDir}/bin"
params.launch_dir = "${workflow.launchDir}"
params.tru_seq_pe_fasta_path = "${workflow.launchDir}/TruSeq3-PE.fa"

// Originally we were downloading the fastqs but this was too unrelible.
// There were constant time outs and other errors. As such I think it is
// safest to start from a folder with fastq.gz file in it
if (params.from_download){
    Channel.fromList(params.sra_list).set{ch_download_fastq}

    // We are getting an error when trying to initialize all 10 of the fastq-dump requests at once.
    // when running in tmux. But this doesn't seem to happen outside of tmux
    process download_fastq{
        tag "${sra_id}"
        
        publishDir path: "raw_reads", mode: "copy"

        input:
        val sra_id from ch_download_fastq

        output:
        tuple file("*_1.fastq.gz"), file("*_2.fastq.gz") into ch_correct_fastq_name_format_input

        script:
        """
        python3 ${params.bin_dir}/ena-fast-download.py $sra_id
        """
    }

    // Because we are now using this ena-fast-download
    // we don't get to format the def line of the sequencing files
    // in the same way that we coud with fastq-dump.
    // The format required is quite specific and trimmomatic 
    // will throw us errors if it is not correct.
    // Here we will remove the SRRR part of the name before the space
    // And re-write out the file.
    process correct_fastq_name_format{
        tag "$fastq_one"
        input:
        tuple file(fastq_one), file(fastq_two) from ch_correct_fastq_name_format_input

        output:
        tuple file("*_nf_1.fastq.gz"), file("*_nf_2.fastq.gz") into ch_rename_fastqs

        script:
        """
        gunzip -f $fastq_one
        gunzip -f $fastq_two
        python3 ${params.bin_dir}/correct_fastq_name_format.py
        gzip *.fastq
        """
    }

    // Unfortunately in the above process we can't output the files with exaclty the same name as they
    // are not detected by the search of output names.
    // Becuase we want them to have the same name as the files
    // that might have already been in a downloaded directory
    // we will quickly modify the name in this process
    process rename_fastq_files{
        tag "$fastq_one"

        input:
        tuple file(fastq_one), file(fastq_two) from ch_rename_fastqs
        
        output:
        tuple file("*_1.fastq.gz"), file("*_2.fastq.gz") into ch_fastqc_pre_trim, ch_trimmomatic_input

        script:
        new_name_one = fastq_one.getName().replaceAll('_nf_1.fastq.gz' ,'_1.fastq.gz')
        new_name_two = fastq_two.getName().replaceAll('_nf_2.fastq.gz' ,'_2.fastq.gz')
        """
        mv $fastq_one $new_name_one
        mv $fastq_two $new_name_two
        """
    }

}else{
    // Start from a directory that already contains fastq.gz files in paried sets
    // This will pair the SRR base with the tuple of the files.
    // To match the format of the ch_make_trimmomatic_input
    // We need to process this so that we get rid of the SRR component
    // Probably best if we use map for this.
    Channel.fromFilePairs("${params.raw_reads_dir}/SRR*_{1,2}.fastq.gz").map{it[1]}.into{ch_fastqc_pre_trim; ch_trimmomatic_input}
}



// // Now fastqc the files
// // This does work
// process fastqc_pre_trim{
//     tag "${fastq_file}"
//     conda "envs/nf_general.yaml"
//     storeDir "nf_fastqc_pre_trim"

//     input:
//     file fastq_file from ch_fastqc_pre_trim.flatten()

//     output:
//     file "*_fastqc.html" into ch_fastqc_pre_trim_output

//     script:
//     """
//     fastqc -o . $fastq_file
//     """
// }

// Now trim the files
// NB wherever possible it is really important to be as exact with the output files
// as possible. For example, here we go to the trouble of being specific about th SRR
// base instead of just *1_.fq.gz. The lack of specificity is OK if we are using the
// '-resume' cache system, because you only look for the output files in each of 
// the directories of the individual process instances. However, when we start to use
// storeDir, which is super useful, to do cacheing, the lack of specificity beomes a problem.
// This is because when we come to read in the output files, we are doing so from the directory
// that holds all outputs from all samples. This read in happens for each instance of the process
// that is being skipped and so for each instance we end up with the contents of the directory.
// By being specific we ensure that only the files that would ahve realted to the particular instance
// i.e. to a praticular sample, are read in.
process trimmomatic{
	tag "${fastq_file_one}"
    conda "envs/nf_general.yaml"
	storeDir "nf_trimmed"

	input:
	tuple file(fastq_file_one), file(fastq_file_two) from ch_trimmomatic_input
	
	output:
	// Output that will be used for the error_correction
	// This is a list of tuples that are the 1P and 2P output files only
	tuple file("${fastq_file_one.getName().replaceAll('_1.fastq.gz', '')}*1P.fq.gz"), file("${fastq_file_one.getName().replaceAll('_1.fastq.gz', '')}*2P.fq.gz") into ch_rcorrect_input

	script:
	outbase = fastq_file_one.getName().replaceAll('_1.fastq.gz', '.trimmed.fq.gz')
	"""
	trimmomatic PE -threads ${params.trimmomatic_threads} -basein ${fastq_file_one} \\
		-baseout $outbase \\
		ILLUMINACLIP:${params.tru_seq_pe_fasta_path}:2:30:10:2:keepBothReads \\
		LEADING:3 TRAILING:3 MINLEN:36 HEADCROP:11
	"""
}

ch_rcorrect_input.view()

// // Now error correction
// process rcorrector{
//     tag "${trimmed_read_one}"
//     conda "envs/nf_general.yaml"
//     cpus params.rcorrector_threads
    
//     storeDir "nf_error_corrected"
    
//     input:
//     tuple file(trimmed_read_one), file(trimmed_read_two) from ch_rcorrect_input

//     output:
//     tuple file("${trimmed_read_one.getName().replaceAll('.trimmed_1P.fq.gz', '')}*1P.cor.fq.gz"), file("${trimmed_read_one.getName().replaceAll('.trimmed_1P.fq.gz', '')}*2P.cor.fq.gz") into ch_trinity_input

//     script:
//     """
//     run_rcorrector.pl -1 $trimmed_read_one -2 $trimmed_read_two -od . -t ${task.cpus}
//     """
// }




// // Now do the trinity assembly
// // NB I have changed the code in here subtly without re running it in order to implement the
// // storeDir directive. I have added a final mv line to the script to rename the ambiguous trinity fasta
// // file so that the name is specific to the sample SRR base name.
// // I have also changed the output dir from being a specific folder for each of the samples (due to the 
// // ambiguity in the previous naming system) to all being held in the nf_trinity_assembly dir.
// // I will make the changes by hand now to the outputs that are already in this directory.
// process trinity{
//     cache 'lenient'
//     tag "${corrected_read_one}"
//     conda "envs/nf_general.yaml"
//     cpus params.trinity_threads
//     storeDir "nf_trinity_assembly"

//     input:
//     tuple file(corrected_read_one), file(corrected_read_two) from ch_trinity_input

//     output:
//     file "${corrected_read_one.getName().replaceAll('.trimmed_1P.cor.fq.gz', '')}*Trinity.fasta" into ch_remove_short_iso_forms_trinity_input
    
//     script:
//     // NB that the output directory for each trinity assembly must have 'trinity' in it.
//     """
//     Trinity --left $corrected_read_one --right $corrected_read_two --seqType fq --max_memory 150G --CPU ${task.cpus} \\
//     --min_contig_length 250 --output trinity --full_cleanup
//     mv trinity.Trinity.fasta ${corrected_read_one.getName().replaceAll('.trimmed_1P.cor.fq.gz','')}.trinity.Trinity.fasta
//     """
// }

// // We will have to write a process for including the other transcriptomes here.
// // Basically, we will do a collect on the ch_remove_short_iso_forms_trinity_input channel
// // And then get merge this with a read in of two file paths
// // process add_new_trinity_assemblies{
// //     cache 'lenient'
// //     tag "add_new_trinity"
// //     conda "envs/nf_general.yaml"

// //     input:

// // }

// // // Now remove the short isos from the trinity assembly
// // // NB we were having problems getting the cache of this process to work.
// // // The problem was that I was using $workflow.launch_dir as the base to the /bin/remove... path
// // // However, when I looked at the hash logs, the value of this ($workflow.launch_dir) was giving a different hash
// // // hash each time. I think this was because what was getting hashed was the state of this directory.
// // // Obviously the state was changing very often.
// // // To get the hash to work in the end I have replaced the $workflow.launch_dir with a params.bin_dir
// // // variable that I have set at the beginning of this file using the $workflow.launch_dir variable.
// // // I will look to see now whether the hash logs are looking at this new variable as a string or at
// // // the state of this directory. If they're looking at the state of the directory then caches may
// // // fail when we add additional script files to the /bin directory.
// // This is where we will need to incorporate the new trinity assemblies.
// // I have tried to add the two additional transcriptomes by collecting mixing and flattening
// process remove_short_isos{
//     cache 'lenient'
//     tag "${trinity_assembly_fasta}"

//     conda "envs/nf_python_scripts.yaml"
//     storeDir "nf_trinity_assembly"

//     input:
//     file trinity_assembly_fasta from ch_remove_short_iso_forms_trinity_input.collect().mix(Channel.fromPath(["/home/humebc/projects/parky/breviolum_transcriptomes/nf_trinity_assembly/BreviolumB5_Sradians.trinity.Trinity.fasta", "/home/humebc/projects/parky/breviolum_transcriptomes/nf_trinity_assembly/Breviolumfaviinorum_Pclivosa.trinity.Trinity.fasta"])).flatten()

//     output:
//     file "${trinity_assembly_fasta.getName().replaceAll('.trinity.Trinity.fasta','')}*.long_iso_only.fasta" into ch_orf_prediction_input
    
//     script:
//     """
//     python3 ${params.bin_dir}/remove_short_isos.py $trinity_assembly_fasta
//     """
// }


// // Do ORF prediction using transdecoder
// // Concurrently rename the default names output by transdecoder
// // To make the long_iso_orf files unique and related to their transcriptome we will append the srrname
// // NB 08/03/2020 I am going to modify the flow a little to reduce the redundancy.
// // I seem to be passing around the .cds and .pep file when only the .pep file is required.
// // The .cds is not required until later.
// process orf_prediction{
//     cache 'lenient'
//     tag "${srrname}"
//     conda "envs/nf_transdecoder.yaml"
//     storeDir "nf_transdecoder"
    
//     input:
//     file long_iso_trinity from ch_orf_prediction_input

//     output:
//     file "${long_iso_trinity.getName().replaceAll('.trinity.Trinity.long_iso_only.fasta','')}*.pep" into ch_remove_multi_orfs_input
//     file "${long_iso_trinity.getName().replaceAll('.trinity.Trinity.long_iso_only.fasta','')}*.cds" into ch_write_unaligned_cds_fastas_fas_input
    
//     script:
//     """
//     TransDecoder.LongOrfs -t $long_iso_trinity -O .
//     mv longest_orfs.cds ${long_iso_trinity.getName().replaceAll('.trinity.Trinity.long_iso_only.fasta','')}_longest_iso_orfs.cds
//     mv longest_orfs.pep ${long_iso_trinity.getName().replaceAll('.trinity.Trinity.long_iso_only.fasta','')}_longest_iso_orfs.pep
//     """
// }


// // The transdecoder output can have multiple ORFs predicted per transcript
// // We will once again only keep one representative per transcript and work with this for the
// // ortholog prediction
// // Sonic Parnoid runs from a single directory containing all of the fastas
// // To enable this we will publish each of the fastas into a single directory
// process remove_multi_orfs_from_pep{
//     cache 'lenient'
//     tag "${pep_file}"
//     conda "envs/nf_python_scripts.yaml"
//     storeDir "nf_sonicparanoid"

//     input:
//     file pep_file from ch_remove_multi_orfs_input

//     output:
//     //tuple file("*.single_orf.pep"), file(cds_file) into ch_sonic_paranoid_input
//     file("${pep_file.getName().replaceAll('_longest_iso_orfs.pep','')}*.single_orf.pep") into ch_sonicparanoid_input
    
//     script:
//     output_path = pep_file.getName().replaceAll("longest_iso_orfs.pep", "longest_iso_orfs.single_orf.pep")
    
//     "python3 ${params.bin_dir}/unique_orfs_from_pep.py $pep_file"
// }


// process sonicparanoid{
//     cache 'lenient'
//     tag "sonicparanoid"
//     cpus params.sonicparanoid_threads
//     conda "envs/nf_sonicparanoid.yaml"
//     storeDir "nf_sonicparanoid/output_12"
    
//     input:
//     // We won't actually use this input. It is just here to link
//     // The processes
//     file pep_file from ch_sonicparanoid_input.collect()

//     output:
//     file "single-copy_groups.tsv" into ch_screen_sonicparanoid_output_input

//     script:
//     """
//     sonicparanoid -i ${params.launch_dir}/nf_sonicparanoid -o . -t ${task.cpus}
//     find . -name single-copy_groups.tsv | xargs -I {} mv {} .
//     """
// }


// // The sonic paranoid output table contains orthologs that were not found in all of the transciptomes
// // We will drop these orthologs and write out the .tsv again
// process screen_sonicparnoid_output{
//     tag "screen sonicparanoid"
//     conda "envs/nf_python_scripts.yaml"
//     storeDir "nf_sonicparanoid/output_12"

//     input:
//     file single_copy_tsv from ch_screen_sonicparanoid_output_input

//     output:
//     file "screened_orthologs.tsv" into ch_write_unaligned_cds_fastas_tab_input

//     script:
//     """
//     python3 ${params.bin_dir}/screen_orthologs.py $single_copy_tsv
//     """
// }


// // Work through the screened_orthologs and for each ortholog make a directory
// // and write out a fasta that contains the sequence for each of the strains
// // Input to the script will be the .tsv. Each of the .cds fastas will be
// // in the nextflow working directory due to the channel input
// // The column titles of the screened ortholog .tsv are the .pep file
// // names (e.g. SRR1793320_longest_iso_orfs.single_orf.pep)
// // The cds files that we will be getting the sequences from are all prefixed
// // with the SRRXXX (e.g. SRR1793320_longest_iso_orfs.cds). As such we can use
// // this SRRXXX in the python script to map between the table and the .cds files
// process write_unaligned_cds_fastas{
//     tag "write_unaligned_cds_fastas"
//     conda "envs/nf_python_scripts.yaml"
//     storeDir "nf_local_alignments"
    
//     input:
//     file screened_orth_table from ch_write_unaligned_cds_fastas_tab_input
//     file cds_fastas from ch_write_unaligned_cds_fastas_fas_input.collect()

//     output:
//     file "*_unaligned_cds.fasta" into ch_align_using_guidance_input

//     script:
//     """
//     python3 ${params.bin_dir}/write_out_unaligned_cds_fastas.py $screened_orth_table
//     """
// }

// // Align each of the unaligned cds files.
// // We will do this in two processes.
// // First we will perform the guidance analysis
// // Then we will use the three output files to create the aligned fastas
// // That have been cropped and had the appropriate low quality columns dropped
// // As always, Guidance is proving tricky to get to run
// // We have had to explicityly add per-bioperl to the env yaml
// // Then, the bash wrapper for the perl execution of the guidance.pl
// // is not working so you have to explicityly run perl and the full path to guidance.pl
// // Because we can't know what the full path is, we will write a quick python script
// // to find this out and run it.
// // Guidance seems to be quite unstable. It seems to produce errors for no reason.
// // However, the pipeline can simply be restarted and everything seems to behave OK.
// // As such I have now implemented an automatic retry process where the python
// // script will attempt to re-run the guidance analysis in question before quitting.
// process align_using_guidance{
//     tag "${unaligned_fasta.toString().split('_')[0]}"
//     conda "envs/nf_guidance.yaml"
//     storeDir "nf_local_alignments"
//     input:
//     file unaligned_fasta from ch_align_using_guidance_input.flatten()

//     output:
//     tuple file("${unaligned_fasta.toString().split('_')[0]}*.MAFFT.Guidance2_res_pair_res.PROT.scr"), file("${unaligned_fasta.toString().split('_')[0]}*.MAFFT.PROT.aln"), file("${unaligned_fasta.toString().split('_')[0]}*.MAFFT.aln.With_Names") into ch_process_guidance_output_input

//     script:
//     orth_group_id = unaligned_fasta.toString().split('_')[0]
//     """
//     python3 ${params.bin_dir}/run_guidance.py $unaligned_fasta
//     """
// }

// // This process will use the three ouput files from align_using_guidance
// process process_guidance_output{
//     tag "${aa_cols_score_file_path.toString().split("/")[-1].split(/\./)[0]}"
//     storeDir "nf_local_alignments"

//     input:
//     tuple file(aa_cols_score_file_path), file(aa_alignment_file_path), file(cds_alignment_file_path) from ch_process_guidance_output_input

//     output:
//     file "${aa_cols_score_file_path.toString().split("/")[-1].split(/\./)[0]}*_cropped_aligned_aa.fasta" into ch_model_test_input
//     file "${aa_cols_score_file_path.toString().split("/")[-1].split(/\./)[0]}*_cropped_aligned_cds.fasta" into ch_run_codeml_align_input

//     script:
//     """
//     python3 ${params.bin_dir}/process_guidance_output.py $aa_cols_score_file_path $aa_alignment_file_path $cds_alignment_file_path
//     """
// }


// // Now find the best protein evo model
// // In some cases the fasta files that resulted from the process_guidance_output may be empty
// // As such we will wrap the modeltest-ng running in a python script that will check for this.
// // If there is this problem with the aligned fasta then we will return code 0
// // before making sure to delete any *.out file that might have been made
// // We will make it so that the *_prottest_result.out is output empty.
// // We can then check for the empty .out file in the make_master_alignment
// process model_test{
//     tag "${cropped_aligned_aa_fasta.toString().split('_')[0]}"
//     conda "envs/nf_modeltest-ng.yaml"
//     storeDir "nf_prot_out"

//     input:
//     file cropped_aligned_aa_fasta from ch_model_test_input

//     output:
//     tuple file("${cropped_aligned_aa_fasta.toString().split('_')[0]}*_prottest_result.out"), file(cropped_aligned_aa_fasta) into ch_make_master_alignment_input

//     script:
//     """
//     python3 ${params.bin_dir}/run_model_test.py $cropped_aligned_aa_fasta
//     """
// }

// // // ch_make_master_alignment_input.collect().subscribe {  println "Got: $it"  }

// // To make the master tree we will work witha single process that
// // will need to iterthrough each of the protein model outputs.
// // It will also need access to the aa cropped and alignment files
// // We will supply both of these in two seperate input channels
// // The output will be a master fasta and a q file for raxml that delimits the partitions
// // that can then be fed into the treemaking
// process make_master_alignment_and_q_file{
//     tag "make_master_alignment"
//     conda "envs/nf_python_scripts.yaml"
//     storeDir "nf_master_fasta_and_q_file"

//     input:
//     file out_and_aligned_fasta_files from ch_make_master_alignment_input.collect()

//     output:
//     tuple file("master_fasta_for_tree.fasta"), file("q_partition_file.q") into ch_make_tree_input

//     script:
//     """
//     python3 ${params.bin_dir}/make_master_alignment.py ${params.launch_dir}/nf_prot_out
//     """
// }

// // Make a ML tree using raxml
// process make_tree{
//     tag "make_tree"
//     cpus params.raxml_threads
//     conda "envs/nf_raxml.yaml"
//     storeDir "nf_master_tree"

//     input:
//     tuple file(master_fasta), file(q_partition_file) from ch_make_tree_input

//     output:
//     tuple file("RAxML_bestTree.strain_dn_ds"), file("RAxML_bipartitionsBranchLabels.strain_dn_ds"), file("RAxML_bipartitions.strain_dn_ds") into ch_annotate_tree_input
//     file "RAxML_bestTree.strain_dn_ds" into ch_run_codeml_tree_input
//     script:
//     """
//     raxmlHPC-PTHREADS-AVX2 -s $master_fasta -q $q_partition_file -x 183746 -f a, -p \\
//     83746273 -# 1000 -T ${task.cpus} -n strain_dn_ds -m PROTGAMMAWAG
//     """
// }

// // Apply more human readable labels to the tree
// // Currently the labels are the SRRXXX this will
// //add the strain and species to the labels
// process annotate_tree{
//     tag "annotate_tree"
//     conda "envs/nf_python_scripts.yaml"
//     storeDir "nf_master_tree"

//     input:
//     tuple file(tree_one), file(tree_two), file(tree_three) from ch_annotate_tree_input

//     output:
//     tuple file("RAxML_bestTree.strain_dn_ds_named"), file("RAxML_bipartitionsBranchLabels.strain_dn_ds_named"), file("RAxML_bipartitions.strain_dn_ds_named") into ch_tree_for_dnds_output

//     script:
//     """
//     python3 ${params.bin_dir}/annotate_tree.py $tree_one $tree_two $tree_three
//     """
// }

// // Here we create codeml control files.
// // Previously we were going to a lot of effort to incorporate
// // some sort of parallelisation. However, let's see if we can make use of nextflow
// // format here and run one instance per cds input.
// // When we run this we will check to see that there are sequences in the alignment and 
// // that the alignment is divisible by 3. If either of these assertions fails
// // we will exit without error. This means that the *.out file output needs to be optional
// // Having an optional *.out file causes issues for storeDir as it causes all process to skip
// // rather than look for a .out file. To work around this we have written in the output
// // of a file called status.txt that will hold a '0' or a '1' depending on whether the
// // process caused an error or not. This way the process will be forced to run
// process run_codeml{
//     tag "$cds_file"
//     conda "envs/nf_codeml.yaml"
//     storeDir "nf_codeml_out"
//     input:
//     tuple file(cds_file), file(tree) from ch_run_codeml_align_input.combine(ch_run_codeml_tree_input)
    
//     output:
//     file "${cds_file.toString().split('_')[0]}*_codeml_results.out" optional true into ch_collate_codeml_results_intput
//     file "${cds_file.toString().split('_')[0]}_status.txt" into ch_force_process_run_output
//     script:
//     """
//     python3 ${params.bin_dir}/run_codeml.py $cds_file $tree
//     """
// }

// process collate_codeml_results{
//     tag "collate_codeml_results"
//     conda "envs/nf_python_scripts.yaml"
//     storeDir "nf_dnds_df_summary"

//     input:
//     file codeml_out_file from ch_collate_codeml_results_intput.collect()
//     output:
//     file "codeml_results_df.csv" into ch_collate_codeml_results

//     script:
//     """
//     python3 ${params.bin_dir}/collate_codeml_results.py ${params.sra_list_as_csv}
//     """
// }