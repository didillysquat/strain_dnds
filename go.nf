#!/usr/bin/env nextflow
// This nextflow doc will be concerned with performing the GO anotations
// We want to annotate every predicted ORF that was used as input to the
// sonicparanoid input.
// We will try working on a per transciptome basis
// As such there will be 14 transcriptomes we will want to work with
// The individual 12 in the 10_2 analysis and he additional 2 that are the
// b_minitum and b_pseudo combined transciptomes used in the '4' analyses.
// We were originally doing this work with the blast algorithm
// but this was taking an enormous amount of time. As such we are going to 
// try switching to the mmseqs algo. As a comparison we will check the outputs
// of the swissprot blast done with blast against the swissprot blast done with
// mmseqs 2, just as a sanity check.
params.pep_input_directories = "/home/humebc/projects/parky/strain_dnds/nf_sonicparanoid/10_2/,/home/humebc/projects/parky/strain_dnds/nf_sonicparanoid/4/"
params.path_to_uniprot_sprot_mmseqs = "/share/databases/uniprot_sprot/mmseqs_uniprot_sprot/mmseqs_uniprot_sprot.targetDB"
params.path_to_uniprot_trembl_mmseqs = "/share/databases/uniprot_trembl/mmseqs_uniprot_trembl/mmseqs_uniprot_trembl.targetDB"
params.path_to_nr_mmseqs = "/share/databases/nr/mmseqs_nr/nr.targetDB"
params.nf_mmseqs_trembl_query_dbs = "/home/humebc/projects/parky/strain_dnds/nf_mmseqs_trembl_query_dbs"
params.nf_mmseqs_nr_query_dbs = "/home/humebc/projects/parky/strain_dnds/nf_mmseqs_nr_query_dbs"
Channel.fromPath(["/home/humebc/projects/parky/strain_dnds/nf_sonicparanoid/10_2/*.pep", "/home/humebc/projects/parky/strain_dnds/nf_sonicparanoid/4/SRR_b_min_c_longest_iso_orfs.single_orf.pep", "/home/humebc/projects/parky/strain_dnds/nf_sonicparanoid/4/SRR_b_psyg_c_longest_iso_orfs.single_orf.pep"]).into{ch_input_swiss_prot_blast; ch_input_make_mmseqs_query_dbs_sprot}
params.bin_dir = "${workflow.launchDir}/bin"
params.goa_uniprot_all_gaf_path = "/share/databases/goa_uniprot_all/goa_uniprot_all.gaf"
params.go_cache_dir = "/home/humebc/projects/parky/strain_dnds/nf_go_annotations/cache"

// Do the BLAST swiss_prot_alignments
process swiss_prot_blast{
    cache 'lenient'
    tag "${pep_file.toString().replaceAll('.pep', '')}"
    storeDir 'nf_blast_results'
    conda "envs/nf_go.yaml"
    cpus params.blast_sprot_threads

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

// The MMseqs algo takes a query db, a target db and outputs a results db.
// The target DBs, we will make outside of this work flow and they will reside in the
// zygote databases dir. For the input dbs, we will want to make one for each of the 14 transcriptomes we are working with

// base_name = "${pep_file.getName().replaceAll('_longest_iso_orfs.single_orf.pep', '.queryDB')}"
//     tuple file($base_name), \
    // file("${pep_file.toString().replaceAll('_longest_iso_orfs.single_orf.pep', '.queryDB')}.dbtype"), \
    // file("${pep_file.toString().replaceAll('_longest_iso_orfs.single_orf.pep', '.queryDB')}_h"), \
    // file("${pep_file.toString().replaceAll('_longest_iso_orfs.single_orf.pep', '.queryDB')}_h.dbtype"), \
    // file("${pep_file.toString().replaceAll('_longest_iso_orfs.single_orf.pep', '.queryDB')}_h.index"),  \
    // file("${pep_file.toString().replaceAll('_longest_iso_orfs.single_orf.pep', '.queryDB')}.index"), \
    // file("${pep_file.toString().replaceAll('_longest_iso_orfs.single_orf.pep', '.queryDB')}.lookup"), \
    // file("${pep_file.toString().replaceAll('_longest_iso_orfs.single_orf.pep', '.queryDB')}.source"), into ch_input_mmseqs_swissprot


process make_mmseqs_query_dbs_sprot{
    cache 'lenient'
    tag "${pep_file.toString().replaceAll('.pep', '')}"
    storeDir 'nf_mmseqs_sprot_query_dbs'
    conda "envs/nf_go.yaml"

    input:
    file pep_file from ch_input_make_mmseqs_query_dbs_sprot

    // We will put all of the output file into the publishDir
    output:
    tuple file("${pep_file.toString().replaceAll('_longest_iso_orfs.single_orf.pep', '.queryDB')}"),\
    file("${pep_file.toString().replaceAll('_longest_iso_orfs.single_orf.pep', '.queryDB')}.dbtype"),\
    file("${pep_file.toString().replaceAll('_longest_iso_orfs.single_orf.pep', '.queryDB')}_h"),\
    file("${pep_file.toString().replaceAll('_longest_iso_orfs.single_orf.pep', '.queryDB')}_h.dbtype"),\
    file("${pep_file.toString().replaceAll('_longest_iso_orfs.single_orf.pep', '.queryDB')}_h.index"),\
    file("${pep_file.toString().replaceAll('_longest_iso_orfs.single_orf.pep', '.queryDB')}.index"),\
    file("${pep_file.toString().replaceAll('_longest_iso_orfs.single_orf.pep', '.queryDB')}.lookup"),\
    file("${pep_file.toString().replaceAll('_longest_iso_orfs.single_orf.pep', '.queryDB')}.source") into ch_input_mmseqs_swissprot

    script:
    out_name = "${pep_file.getName().replaceAll('_longest_iso_orfs.single_orf.pep', '.queryDB')}"
    """
    mmseqs createdb $pep_file "${pep_file.getName().replaceAll('_longest_iso_orfs.single_orf.pep', '.queryDB')}" --createdb-mode 1
    """
}

// ch_input_mmseqs_swissprot.view()

process mmseqs_search_swissprot{
    cache 'lenient'
    cpus params.mmseqs_sprot_threads
    tag "$base_db_file"
    storeDir 'nf_mmseqs_sprot_results'
    conda "envs/nf_go.yaml"
    
    input:
    tuple file(base_db_file), file(dbtype), file(h_file), file(h_dbtype), file(h_index), file(index_file), file(lookup_file), file(source_file) from ch_input_mmseqs_swissprot

    output:
    file "${base_db_file.getName().replaceAll('.queryDB', '.mmseqs.out.blast_format')}" into ch_input_mmseqs_trembl_in_pep_dbs

    script:
    resultsdb_name = "${base_db_file.getName().replaceAll('.queryDB', '.resultdb')}"
    blastformat_result_name = "${base_db_file.getName().replaceAll('.queryDB', '.mmseqs.out.blast_format')}"
    """
    mmseqs search $base_db_file ${params.path_to_uniprot_sprot_mmseqs} $resultsdb_name tmp --threads ${task.cpus} -e 0.0000000001 -v 3
    mmseqs convertalis $base_db_file ${params.path_to_uniprot_sprot_mmseqs} $resultsdb_name $blastformat_result_name --threads ${task.cpus} --format-output query,target,evalue,pident,qcov -v 3
    """
}


// We have generated a python script to check the agreement between the two outputs
// It is called mmseqs_blast_agreement. It checks how many of the blast results
// are in agreement with the mmseqs results for the swissport querying.
// It shows us that approximately 70% of the returned mmseqs results are in perfect agreement
// (i.e. fist best match) with the first best match of the BLAST results.
// If considering the first 5 blast results, then you get agreement up to 96%.
// Likely, agreement could be increased if you considered the top 5 of the blast
// and the top 5 of the mmseqs, but we did not do this. I think the 70% result is 
// easily good enough to rely on mmseqs instead of BLAST. We also tested whether increasing
// from the default sensitivity value of mmseqs made a big difference (-s). It did not,
// so we will stick with the default value.


// For each of the out files, we'll want to look to see which of the sequeneces got a search hit
// that we can associate with a GO annotation. If we cannot associate to a GO annotation
// then we will need to move this sequence into the new .pep to go into trembl
process mmseqs_write_out_new_in_pep_for_trembl_search{
    cache 'lenient'
    conda "envs/nf_go.yaml"
    tag "${mmseqs_out.toString().replaceAll('.mmseqs.out.blast_format', '')}"
    storeDir params.nf_mmseqs_trembl_query_dbs
    
    input:
    file mmseqs_out from ch_input_mmseqs_trembl_in_pep_dbs

    output:
    file "$new_pep_file" into ch_input_mmseqs_trembl_make_in_dbs
    
    script:
    // We will write out the new .pep file to the local work directory, and then run the mmseqs createdb
    // function here and pull out the new db files.
    base_name = "${mmseqs_out.toString().replaceAll('.mmseqs.out.blast_format', '')}"
    pep_file = "${base_name}_longest_iso_orfs.single_orf.pep"
    new_pep_file = "${base_name}.mmseqs.trembl.in.pep"
    """
    python3 ${params.bin_dir}/write_out_trembl_in_pep.py $mmseqs_out "${params.pep_input_directories}" $pep_file ${params.goa_uniprot_all_gaf_path} ${params.go_cache_dir}
    """
}

process mmseqs_trembl_make_in_dbs{
    cache 'lenient'
    conda "envs/nf_go.yaml"
    tag "${pep_file.toString().replaceAll('.mmseqs.out.blast_format', '')}"
    storeDir params.nf_mmseqs_trembl_query_dbs
    
    input:
    file pep_file from ch_input_mmseqs_trembl_make_in_dbs

    output:
    tuple file("${query_db_name}"),\
    file("${query_db_name}.dbtype"),\
    file("${query_db_name}_h"),\
    file("${query_db_name}_h.dbtype"),\
    file("${query_db_name}_h.index"),\
    file("${query_db_name}.index"),\
    file("${query_db_name}.lookup"),\
    file("${query_db_name}.source") into ch_input_mmseqs_trembl_search
    
    script:
    // We will write out the new .pep file to the local work directory, and then run the mmseqs createdb
    // function here and pull out the new db files.
    query_db_name = "${pep_file.toString().replaceAll('.pep', '.querydb')}"
    """
    mmseqs createdb $pep_file $query_db_name --createdb-mode 1
    """
}

process mmseqs_search_trembl{
    cache 'lenient'
    cpus params.mmseqs_trembl_threads
    tag "$base_db_file"
    storeDir 'nf_mmseqs_trembl_results'
    conda "envs/nf_go.yaml"
    
    input:
    tuple file(base_db_file), file(dbtype), file(h_file), file(h_dbtype), file(h_index), file(index_file), file(lookup_file), file(source_file) from ch_input_mmseqs_trembl_search

    output:
    file "$blastformat_result_name" into ch_output_mmseqs_search_trembl

    script:
    resultsdb_name = "${base_db_file.getName().replaceAll('.querydb', '.resultdb')}"
    blastformat_result_name = "${base_db_file.getName().replaceAll('.in.querydb', '.out.blast_format')}"
    """
    mmseqs search $base_db_file ${params.path_to_uniprot_trembl_mmseqs} $resultsdb_name tmp --threads ${task.cpus} -e 0.0000000001 -v 3
    mmseqs convertalis $base_db_file ${params.path_to_uniprot_trembl_mmseqs} $resultsdb_name $blastformat_result_name --threads ${task.cpus} --format-output query,target,evalue,pident,qcov -v 3
    """
}

// // Once the trembl search is complete, we will want to look through the output files to get an idea
// // of how many of the query seqs still haven't had a hit in either the sprot or trembl databases.
// // We can then work with these sequences to either put them through the nr database or just 
// // link them up to GO annotations.
// // In Seb's previous paper, he ran the sequeneces against nr, but I don't see how to get from
// // the nr matches to a list of GO terms. Unless we use something like interproscan.

process mmseqs_write_out_new_in_pep_for_nr_search{
    cache 'lenient'
    conda "envs/nf_go.yaml"
    tag "${mmseqs_out.toString().replaceAll('.mmseqs.trembl.out.blast_format', '')}"
    storeDir params.nf_mmseqs_nr_query_dbs
    
    input:
    file mmseqs_out from ch_output_mmseqs_search_trembl

    output:
    file "$new_pep_file" into ch_output_mmseqs_in_nr_pep
    
    script:
    // We will write out the new .pep file to the local work directory, and then run the mmseqs createdb
    // function here and pull out the new db files.
    base_name = "${mmseqs_out.toString().replaceAll('.mmseqs.trembl.out.blast_format', '')}"
    old_pep_file = "${base_name}.mmseqs.trembl.in.pep"
    new_pep_file = "${base_name}.mmseqs.nr.in.pep"
    """
    python3 ${params.bin_dir}/write_out_nr_in_pep.py $mmseqs_out ${params.nf_mmseqs_trembl_query_dbs} $old_pep_file
    """
}

// We now have a means of referring a nr hit to a GO.
// We match the nr hit to a gene_id hit. Using /share/databases/gene2accession/gene2accession
// We then match the gene_id hit to the GO. Using //share/databases/gene2go/gene2go

process mmseqs_nr_make_in_dbs{
    cache 'lenient'
    conda "envs/nf_go.yaml"
    tag "${pep_file.toString().replaceAll('.mmseqs.nr.in.pep', '')}"
    storeDir params.nf_mmseqs_nr_query_dbs
    
    input:
    file pep_file from ch_output_mmseqs_in_nr_pep

    output:
    tuple file("${query_db_name}"),\
    file("${query_db_name}.dbtype"),\
    file("${query_db_name}_h"),\
    file("${query_db_name}_h.dbtype"),\
    file("${query_db_name}_h.index"),\
    file("${query_db_name}.index"),\
    file("${query_db_name}.lookup"),\
    file("${query_db_name}.source") into ch_input_mmseqs_nr_search
    
    script:
    // We will write out the new .pep file to the local work directory, and then run the mmseqs createdb
    // function here and pull out the new db files.
    query_db_name = "${pep_file.toString().replaceAll('.pep', '.querydb')}"
    """
    mmseqs createdb $pep_file $query_db_name --createdb-mode 1
    """
}

process mmseqs_search_nr{
    cache 'lenient'
    cpus params.mmseqs_trembl_threads
    tag "$base_db_file"
    storeDir 'nf_mmseqs_nr_results'
    conda "envs/nf_go.yaml"
    
    input:
    tuple file(base_db_file), file(dbtype), file(h_file), file(h_dbtype), file(h_index), file(index_file), file(lookup_file), file(source_file) from ch_input_mmseqs_nr_search

    output:
    file "$blastformat_result_name" into ch_output_mmseqs_search_nr

    script:
    resultsdb_name = "${base_db_file.getName().replaceAll('.querydb', '.resultdb')}"
    blastformat_result_name = "${base_db_file.getName().replaceAll('.in.querydb', '.out.blast_format')}"
    """
    mmseqs search $base_db_file ${params.path_to_uniprot_trembl_mmseqs} $resultsdb_name tmp --threads ${task.cpus} -e 0.0000000001 -v 3
    mmseqs convertalis $base_db_file ${params.path_to_uniprot_trembl_mmseqs} $resultsdb_name $blastformat_result_name --threads ${task.cpus} --format-output query,target,evalue,pident,qcov -v 3
    """
}


// Let's assume that we're not going to going further and doing nr searching
// In stead we should read in the rn.in.pep created above and use this as the number of sequences
// that we didn't find matches for. We want to end up with gene ontology annotations for each
// of the uniprot accessions that we found an association with. To do this we will need to
// read in the sprot and trembl mmseqs output files and work with these.
// process get_GO_annotations{
//     cache 'lenient'
//     conda "envs/python_scripts.yaml"
//     tag $base_name
//     storeDir 'nf_go_associations'

//     input:
//     file nr_in_pep from ch_output_mmseqs_in_nr_pep

//     output:
//     file "$annotation_out" into ch_annotation_out

//     script:
//     base_name = "${nr_in_pep.toString.replaceAll('.mmseqs.nr.in.pep', '')}"
//     annotation_out = "${base_name}.go.annotations.csv"
    
//     """
//     python3 ${params.bin_dir}/collect_go_annotations.py $base_name  ${params.pep_input_directories} ${params.nf_mmseqs_trembl_query_dbs} $old_pep_file
//     """
// }