# Produce fastq files for each of the raw transcript sequncing files
# $ snakemake --cores 24 fastqc/b_minutum/SRR17933{20..23}_{1,2}_fastqc.html fastqc/b_psygmophilum/SRR17933{24..27}_{1,2}_fastqc.html
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
rule assemble:
	input:
		fwd="err_corrected/{species}/{sra}.trimmed_1P.cor.fq.gz",
		rev="err_corrected/{species}/{sra}.trimmed_2P.cor.fq.gz"
	output:
		"apple_pie_{species}_{sra}.txt"
	conda:
		"envs/strain_deg.yaml"
	shell:
		"Trinity --left {input.fwd} --right {input.rev} --seqType fq --max_memory 500G --CPU 24 --min_contig_length 250 --output ./trinity_assembly/{wildcards.species}_trinity"
    
