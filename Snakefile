# Produce fastq files for each of the raw transcript sequncing files
# $ snakemake --cores 24 fastqc/b_minutum/SRR17933{20..23}_{1,2}_fastqc.html fastqc/b_psygmophilum/SRR17933{24..27}_{1,2}_fastqc.html
sra_dict = {"b_minutum":["SRR1793320", "SRR1793321", "SRR1793322", "SRR1793323"], "b_psygmophilum": ["SRR1793324", "SRR1793325", "SRR1793326", "SRR1793327"]}
wildcard_constraints:
    sra="SRR\d+"
rule make_fastqc_min:
	input:
		"raw_reads/b_minutum/{sra}_{direction_num}.fastq.gz"
	output:
		"fastqc/b_minutum/{sra}_{direction_num}_fastqc.html"
	shell:
		"fastqc -o fastqc/b_minutum {input}"
rule make_fastqc_psy:
	input:
		"raw_reads/b_psygmophilum/{sra}_{direction_num}.fastq.gz"
	output:
		"fastqc/b_psygmophilum/{sra}_{direction_num}_fastqc.html"
	shell:
		"fastqc -o fastqc/b_psygmophilum {input}"

# Use trimmomatic to remove 3' adapters and and to remove first few bp from begining
# When examining the fastqc the first 6 or so bp looked bad
# snakemake -np trimmed/b_psygmophilum/SRR17933{24..27}.trimmed_1P.fq.gz trimmed/b_minutum/SRR17933{20..23}.trimmed_1P.fq.gz
rule trim_min:
	input:
		"raw_reads/b_minutum/{sra}_1.fastq.gz"
	output:
		"trimmed/b_minutum/{sra}.trimmed_1P.fq.gz",
		"trimmed/b_minutum/{sra}.trimmed_2P.fq.gz",
		"trimmed/b_minutum/{sra}.trimmed_1U.fq.gz",
		"trimmed/b_minutum/{sra}.trimmed_2U.fq.gz"
	threads:6
	shell:
		"trimmomatic PE -threads {threads} -basein {input} -baseout trimmed/b_minutum/{wildcards.sra}.trimmed.fq.gz "
		"ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36 HEADCROP:11"	

rule trim_psy:
	input:
		"raw_reads/b_psygmophilum/{sra}_1.fastq.gz"
	output:
		"trimmed/b_psygmophilum/{sra}.trimmed_1P.fq.gz",
		"trimmed/b_psygmophilum/{sra}.trimmed_2P.fq.gz",
		"trimmed/b_psygmophilum/{sra}.trimmed_1U.fq.gz",
		"trimmed/b_psygmophilum/{sra}.trimmed_2U.fq.gz"
	threads:6
	shell:
		"trimmomatic PE -threads {threads} -basein {input} -baseout trimmed/b_psygmophilum/{wildcards.sra}.trimmed.fq.gz "
		"ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36 HEADCROP:11"	

rule err_correct_min:
	input:
		fwd = "trimmed/b_minutum/{sra}.trimmed_1P.fq.gz",
		rev = "trimmed/b_minutum/{sra}.trimmed_2P.fq.gz"
	output:
		"err_corrected/b_minutum/{sra}.trimmed_1P.cor.fq.gz",
		"err_corrected/b_minutum/{sra}.trimmed_2P.cor.fq.gz"
	threads: 6
	shell:
		"run_rcorrector.pl -1 {input.fwd} -2 {input.rev} -od err_corrected/b_minutum/ -t {threads}"

rule err_correct_psy:
	input:
		fwd = "trimmed/b_psygmophilum/{sra}.trimmed_1P.fq.gz",
		rev = "trimmed/b_psygmophilum/{sra}.trimmed_2P.fq.gz"
	output:
		"err_corrected/b_psygmophilum/{sra}.trimmed_1P.cor.fq.gz",
		"err_corrected/b_psygmophilum/{sra}.trimmed_2P.cor.fq.gz"
	threads: 6
	shell:
		"run_rcorrector.pl -1 {input.fwd} -2 {input.rev} -od err_corrected/b_minutum/ -t {threads}"


