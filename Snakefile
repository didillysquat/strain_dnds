# Produce fastq files for each of the raw transcript sequncing files
# $ snakemake --cores 24 fastqc/b_minutum/SRR17933{20..23}_{1,2}_fastqc.html fastqc/b_psygmophilum/SRR17933{24..27}_{1,2}_fastqc.html
from itertools import combinations
import subprocess
import os
sra_dict = {"b_minutum":["SRR1793320", "SRR1793321", "SRR1793322", "SRR1793323"], "b_psygmophilum": ["SRR1793324", "SRR1793325", "SRR1793326", "SRR1793327"]}
wildcard_constraints:
    sra="SRR\d+"
# Get the fastq files from NCBI using fastq-dump
rule download_sra:
	output:
		"raw_reads/{species}/{sra}_1.fastq",
		"raw_reads/{species}/{sra}_2.fastq"
	conda:
		"envs/strain_deg.yaml"
	shell:
		"fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files {wildcards.sra} -O raw_reads/{wildcards.species}"

# gzip the fastq files downloaded from NCBI using gzip
rule gzip_fastq:
	input:
		"raw_reads/{species}/{sra}_1.fastq",
		"raw_reads/{species}/{sra}_2.fastq"
	output:
		"raw_reads/{species}/{sra}_1.fastq.gz",
		"raw_reads/{species}/{sra}_2.fastq.gz"
	conda:
		"envs/strain_deg.yaml"
	shell:
		"parallel gzip ::: {input}"

rule make_fastqc:
	input:
		"raw_reads/{species}/{sra}_{direction_num}.fastq.gz"
	output:
		"fastqc/{species}/{sra}_{direction_num}_fastqc.html"
	conda:
		"envs/strain_deg.yaml"
	shell:
		"fastqc -o fastqc/{wildcards.species} {input}"

# Use trimmomatic to remove 3' adapters and and to remove first few bp from begining
# When examining the fastqc the first 6 or so bp looked bad
# snakemake -np trimmed/b_psygmophilum/SRR17933{24..27}.trimmed_1P.fq.gz trimmed/b_minutum/SRR17933{20..23}.trimmed_1P.fq.gz
rule trim:
	input:
		"raw_reads/{species}/{sra}_1.fastq.gz"
	output:
		"trimmed/{species}/{sra}.trimmed_1P.fq.gz",
		"trimmed/{species}/{sra}.trimmed_2P.fq.gz",
		"trimmed/{species}/{sra}.trimmed_1U.fq.gz",
		"trimmed/{species}/{sra}.trimmed_2U.fq.gz"
	conda:
		"envs/strain_deg.yaml"
	threads:6
	shell:
		"trimmomatic PE -threads {threads} -basein {input} -baseout trimmed/{wildcards.species}/{wildcards.sra}.trimmed.fq.gz "
		"ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36 HEADCROP:11"
# Use rcorrect to do kmer based error correction of the sequences
# snakemake -np --cores 24 err_corrected/b_minutum/SRR17933{20..23}.trimmed_{1,2}P.cor.fq.gz err_corrected/b_psygmophilum/SRR17933{24..27}.trimmed_{1,2}P.cor.fq.gz
rule err_correct:
	input:
		fwd = "trimmed/{species}/{sra}.trimmed_1P.fq.gz",
		rev = "trimmed/{species}/{sra}.trimmed_2P.fq.gz"
	output:
		"err_corrected/{species}/{sra}.trimmed_1P.cor.fq.gz",
		"err_corrected/{species}/{sra}.trimmed_2P.cor.fq.gz"
	conda:
		"envs/strain_deg.yaml"
	threads: 6
	shell:
		"run_rcorrector.pl -1 {input.fwd} -2 {input.rev} -od err_corrected/{wildcards.species}/ -t {threads}"

# Assemble the QC reads using Trinity
# snakemake -np --cores 24 trinity_assembly/b_minutum/SRR17933{20..23}_trinity/Trinity.fasta trinity_assembly/b_psygmophilum/SRR17933{24..27}_trinity/Trinity.fasta
rule assemble:
	input:
		fwd="err_corrected/{species}/{sra}.trimmed_1P.cor.fq.gz",
		rev="err_corrected/{species}/{sra}.trimmed_2P.cor.fq.gz"
	output:
		"trinity_assembly/{species}/{sra}_trinity.Trinity.fasta",
		"trinity_assembly/{species}/{sra}_trinity.Trinity.fasta.gene_trans_map"
	conda:
		"envs/strain_deg.yaml"
	threads:6
	shell:
		"Trinity --left {input.fwd} --right {input.rev} --seqType fq --max_memory 100G --CPU {threads} "
		"--min_contig_length 250 --output trinity_assembly/{wildcards.species}/{wildcards.sra}_trinity "
		"--full_cleanup"
    

# This rule will create a new trinity fasta file that retains only the longest version of those
# genes that had multiple predicted isoforms
rule remove_short_iso_forms:
	input:
		"trinity_assembly/{species}/{sra}_trinity.Trinity.fasta"
	output:
		"trinity_assembly/{species}/{sra}_trinity.Trinity.long_iso_only.fasta"
	conda:
		"envs/strain_deg.yaml"
	shell:
		"python3.6 scripts/remove_short_isos.py {input} {output}"


# We will use transdecoder to predict ORFs
rule orf_prediction:
	input:
		"trinity_assembly/{species}/{sra}_trinity.Trinity.long_iso_only.fasta"
	output:
		"orf_prediction/{species}/{sra}/longest_iso_orfs.pep"
	shell:
		"TransDecoder.LongOrfs -t {input} -O orf_prediction/{wildcards.species}/{wildcards.sra}"

# The transdecoder output can have multiple ORFs predicted per transcript
# We will once again only keep one representative per transcript and work with this for the
# ortholog prediction
rule remove_multi_orfs_from_pep:
    input:
        "orf_prediction/{species}/{sra}/longest_iso_orfs.pep"
    output:
        "orf_prediction/{species}/{sra}/longest_iso_orfs.single_orf.pep"
    conda:
        "envs/python_scripts.yaml"
    shell:
        "python3.7 scripts/unique_orfs_from_pep.py {input} {output}"


def copy_pep_files():
    # create the sonic paranoid directory if it does not already exist
    dir_to_make = os.path.abspath('sonicparanoid')
    print(f'making dir {dir_to_make}')
    os.makedirs(dir_to_make, exist_ok=True)

    # now copy over each of the files
    for species_key, sra_list in sra_dict.items():
        for sra_val in sra_list:
            from_val = os.path.abspath(f'orf_prediction/{species_key}/{sra_val}/longest_orfs.pep')
            to_val = os.path.abspath(f'sonicparanoid/{sra_val}_longest_orfs.pep')
            print(f'Copying {from_val} to {to_val}')
            subprocess.run(['cp', from_val, to_val])

# sonic paranoid operates on a directory that contains all of the fastafiles that the orthologs will be predicted
# from. So we need to create this directory and copy over the .pep files into this directory
rule copy_pep_files_for_sonicparanoid:
	input:
		expand("orf_prediction/b_minutum/{sra}/longest_iso_orfs.single_orf.pep", sra=sra_dict['b_minutum']),
		expand("orf_prediction/b_psygmophilum/{sra}/longest_iso_orfs.single_orf.pep", sra=sra_dict['b_psygmophilum'])
	output:
		expand("sonicparanoid/{sra}_longest_iso_orfs.single_orf.pep", sra=sra_dict['b_minutum']),
		expand("sonicparanoid/{sra}_longest_iso_orfs.single_orf.pep", sra=sra_dict['b_psygmophilum'])
	run:
		copy_pep_files()

# We will do ortholog prediction using reciprocal blast.
# This was previously done using a combination of inparanoid and multiparanoid.
# This is painstakingly slow and only runs using the old blast.
# Instead we will use sonicparanoid that uses the same algorithms as inparanoid but
# uses mmseq2 instead of blast.
# The conda install of sonicblast doesn't work so you have to install it using pip
# we can work this into the conda env as such
# NB this is acutally causing us problems from within the snakemake file
# we keep getting an ERROR saying that the run ID parkinson_slc was used in a previous run.
"""
channels:
  - bioconda
dependencies:
  - python=3.7
  - pip
  - pip:
      - sonicparanoid==1.2.6
"""
rule orthology_sonic_paranoid:
	input:
		expand("sonicparanoid/{sra}_longest_iso_orfs.single_orf.pep", sra=sra_dict['b_minutum']),
		expand("sonicparanoid/{sra}_longest_iso_orfs.single_orf.pep", sra=sra_dict['b_psygmophilum'])
	output:
		"sonicparanoid/output/runs/parkinson/ortholog_groups/ortholog_groups.tsv"
	log:
		"logs/sonicparanoid.log"
	conda:
		"envs/sonicparanoid.yaml"
	threads:24
	shell:
		"sonicparanoid -i sonicparanoid -o sonicparanoid/output -p parkinson -t {threads}"

rule orthology_sonic_paranoid_slc:
	input:
		expand("sonicparanoid/{sra}_longest_iso_orfs.single_orf.pep", sra=sra_dict['b_minutum']),
		expand("sonicparanoid/{sra}_longest_iso_orfs.single_orf.pep", sra=sra_dict['b_psygmophilum'])
	output:
		"sonicparanoid_out/runs"
	conda:
		"envs/sonicparanoid.yaml"
	threads:24
	shell:
		"sonicparanoid -i sonicparanoid -o sonicparanoid_out -t {threads} -slc"
		# "python3 scripts/install_and_run_sonicparanoid.py"

# Use the output from sonicparanoid to extract those ortholog groups that contain all 8 of the strains
# Also make sure that we only have one ortholog per transcript.
# The output of this will be in the same format as the input
rule extract_unique_cross_strain_orthologs:
	input:
		input_one="sonicparanoid/output/runs/parkinson/ortholog_groups/ortholog_groups.tsv",
		input_two="sonicparanoid/output/runs/parkinson/ortholog_groups/single-copy_groups.tsv"
	output:
		"screened_orthologs/screened_orthologs.tsv"
	shell:
		"python3.6 scripts/screen_orthologs.py {input.input_one} {input.input_two} {output}"

rule tester:
	output:
		  "tester.txt"
	conda:
		 "envs/inparanoid.yaml"
	shell:
		 "sonicparanoid -h"


# ANNOTATION
rule get_swiss_prot_db:
	output:
		"db/swiss_prot/uniprot_sprot.fasta.gz"
	shell:
		"wget -O {output} ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"

rule get_trembl_db:
	output:
		"db/trembl/uniprot_trembl.fasta.gz"
	shell:
		"wget -O {output} ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz"

rule get_ncbi_nr_db:
	output:
		"nr.fasta"
	shell:
		"bash breviolum_transcriptomes/db/ncbi_nr/update_blastdb_wrapper.sh"

rule get_uni_prot_goa_db:
	output:
		"db/uniprot_goa/goa_uniprot_all.gaf.gz"
	shell:
		"wget -O {output} ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz"

rule gunzip_uni_prot_goa_db:
	input:
		"db/uniprot_goa/goa_uniprot_all.gaf.gz"
	output:
		"db/uniprot_goa/goa_uniprot_all.gaf"
	shell:
		"gunzip {input}"

rule gunzip_swiss_prot_db:
	input:
		"db/swiss_prot/uniprot_sprot.fasta.gz"
	output:
		"db/swiss_prot/uniprot_sprot.fasta"
	shell:
		"gunzip {input}"

rule gunzip_tremble_db:
	input:
		"db/trembl/uniprot_trembl.fasta.gz"
	output:
		"db/trembl/uniprot_trembl.fasta"
	shell:
		"gunzip {input}"

rule make_swiss_prot_db:
	input:
		"db/swiss_prot/uniprot_sprot.fasta"
	output:
		"db/swiss_prot/uniprot_sprot.phr",
		"db/swiss_prot/uniprot_sprot.pin",
		"db/swiss_prot/uniprot_sprot.psq"
	shell:
		"makeblastdb -in {input} -out db/swiss_prot/uniprot_sprot -dbtype prot"

rule make_tremble_prot_db:
	input:
		"db/trembl/uniprot_trembl.fasta"
	output:
		"db/trembl/uniprot_trembl.phr",
		"db/trembl/uniprot_trembl.pin",
		"db/trembl/uniprot_trembl.psq"
	shell:
		"makeblastdb -in {input} -out db/trembl/uniprot_trembl -dbtype prot"

# Blast each of the .long_iso_only.fasta files against the swissprot, trembl and ncbi nr databases in that order
# retaining the blast hits that were < 1x10-5. Those that aren't matched at that level blast against the next db.
rule blast_against_swiss_prot:
	input:
		query_fasta = "trinity_assembly/{species}/{sra}_trinity.Trinity.long_iso_only.fasta",
		made_db_path = "db/swiss_prot/uniprot_sprot.psq"
	output:
		"blast_matches/swiss_prot/{species}/{sra}_sp_blast.out.txt"
	threads:6
	shell:
		"blastx -query {input.query_fasta} -out {output} -db db/swiss_prot/uniprot_sprot -num_threads {threads} -outfmt 6 -evalue 1e-5"

# Perform additional blasts if required (i.e. if not all transcripts had a hit. And then create a blast.out.txt that
# is the consolidation of the (upto) 3 blast searches and put in separate directory
#NB this is taking so long that we will abandon doing it for the time being and see if we can
# move forward with the dnds predictions without annotations.
rule additional_blasts:
	# Must ensure that the initial swiss prot blast has been conducted
	# Must also ensure that the other two databases are available
	input:
		swiss_prot_blast_output = "blast_matches/swiss_prot/{species}/{sra}_sp_blast.out.txt",
		long_iso_only_transcript_fasta = "trinity_assembly/{species}/{sra}_trinity.Trinity.long_iso_only.fasta",
		trembl_db_item_path = "db/trembl/uniprot_trembl.pal",
		ncbi_db_item_path = "db/ncbi_nr/nr.pal"
	threads:24 #TODO change this back to 6 when running for real
	output:
		consolidated_blast_out_results = "blast_matches/{species}_consolidated/{sra}_consolidated_blast.out.txt"
	conda:
		"envs/strain_deg.yaml"
	shell:
		"python3.6 scripts/additional_blasts.py {input.swiss_prot_blast_output} {input.long_iso_only_transcript_fasta} "
		"db/trembl/uniprot_trembl db/ncbi_nr/nr "
		"blast_matches/trembl/{wildcards.species}/{wildcards.sra}_trembl_blast.out.txt "
		"blast_matches/ncbi_nr/{wildcards.species}/{wildcards.sra}_ncbi_nr_blast.out.txt "
		"{output.consolidated_blast_out_results} "
		"blast_matches/trembl/{wildcards.species}/{wildcards.sra}_trembl_no_match_blast_in.fasta "
		"blast_matches/ncbi_nr/{wildcards.species}/{wildcards.sra}_ncbi_nr_no_match_blast_in.fasta {threads} "
        "{wildcards.species} {wildcards.sra}"

# We will need a blast database for each of the predicted ORFs.
rule make_db_from_long_iso_only_transcriptomes:
	input:
		"orf_prediction/{species}/{sra}/longest_orfs.pep"
	output:
		"orf_prediction/{species}/{sra}/longest_orfs.phr",
		"orf_prediction/{species}/{sra}/longest_orfs.pin",
		"orf_prediction/{species}/{sra}/longest_orfs.psq"
	shell:
		"makeblastdb -in {input} -dbtype prot -out orf_prediction/{wildcards.species}/{wildcards.sra}/longest_orfs -title longest_orfs"

