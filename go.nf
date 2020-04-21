#!/usr/bin/env nextflow
// This nextflow doc will be concerned with performing the GO anotations
// We want to annotate every predicted ORF that was used as input to the
// sonicparanoid input.
// We will try working on a per transciptome basis
// As such there will be 14 transcriptomes we will want to work with
// The individual 12 in the 10_2 analysis and he additional 2 that are the
// b_minitum and b_pseudo combined transciptomes
// params.pep_input_directories = "/home/humebc/projects/parky/strain_dnds/nf_sonicparanoid/10_2/,/home/humebc/projects/parky/strain_dnds/nf_sonicparanoid/4/"
Channel.fromPath(["/home/humebc/projects/parky/strain_dnds/nf_sonicparanoid/10_2/*.pep", "/home/humebc/projects/parky/strain_dnds/nf_sonicparanoid/4/SRR_b_min_c_longest_iso_orfs.single_orf.pep", "/home/humebc/projects/parky/strain_dnds/nf_sonicparanoid/4/SRR_b_psyg_c_longest_iso_orfs.single_orf.pep"]).set{ch_input_swiss_prot_blast}

process swiss_prot_blast{
    cache 'lenient'
    tag "${pep_file.toString().replaceAll('.pep', '')}"
    publishDir 'nf_blast_results'
    cpus 90

    input:
    file pep_file from ch_input_swiss_prot_blast

    output:
    file "${pep_file.getName().replaceAll('.pep', '')}.sprot.out" into ch_sprot_blast_results_out
    // file "${pep_file.getName().replaceAll('.pep', '')}.trembl_in.fasta" into ch_input_sprot_blast_results_out

    script:
    out_name = "${pep_file.getName().replaceAll('.pep', '')}.sprot.out"
    """
    blastp -query ${pep_file} -db "/share/databases/uniprot_sprot/uniprot_sprot.fasta" -out $out_name \\
    -evalue 1e-10 -num_threads ${task.cpus} -max_target_seqs 10 -outfmt "6 qseqid qacc sseqid sacc evalue pident qcovs"
    """
}

// For each of the out files, we'll want to look to see which of the sequeneces did get a blast hit, and which didn't
// We'll want to put those that didn't into new fasta file that we can then feed into the trembl blast process
// And those that did get a blast hit we'll want to extract the GO annotation from. However, its
// probably easiest if we do the extraction of the GO accession once both the sprot and trembl blasts are completed
// process write_out_trembl_in_fastas{
//     cache 'lenient'
//     tag ""
    
//     input:
//     file blast_out from ch_sprot_blast_results_out

//     output:
//     file "${blast_out.getName().replaceAll('.sprot.out', '.trembl.in.pep')}" into ch_input_trembl_blast

//     script:
//     """
//     python3 ${params.bin_dir}/write_out_trembl_in_fasta.py $blast_out "${params.pep_input_directories}"
//     """
// }