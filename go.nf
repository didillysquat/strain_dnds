#!/usr/bin/env nextflow
// This nextflow doc will be concerned with performing the GO anotations
// We want to annotate every predicted ORF that was used as input to the
// sonicparanoid input.
// We will try working on a per transciptome basis
// As such there will be 14 transcriptomes we will want to work with
// The individual 12 in the 10_2 analysis and he additional 2 that are the
// b_minitum and b_pseudo combined transciptomes

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
    file "${pep_file.getName().replaceAll('.pep', '')}.trembl_in.fasta" into ch_input_sprot_blast_results_out

    script:
    out_name = "${pep_file.getName().replaceAll('.pep', '')}.sprot.out"
    """
    blastp -query ${pep_file} -db "/share/databases/uniprot_sprot/uniprot_sprot.fasta" -out $out_name \\
    -evalue 1e-10 -num_threads ${task.cpus} -max_target_seqs 10 
    """
}