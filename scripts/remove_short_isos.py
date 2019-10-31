"""Read in the trinity fasta that contains all iso forms of given genes.
Remove write out a new fasta that contains only the longest isofrom of each of the genes"""
from Bio import SeqIO
import sys

def remove_short_isoform_sequences():


	input_trinity_fasta_path = sys.argv[1]
	output_trinity_fasta_path = input_trinity_fasta_path.replace('.fasta', '.long_iso_only.fasta')
	print(f'Removing short isoform sequences from: {input_trinity_fasta_path}')
	print('Reading in Trinity assembly')
	trin_file = list(SeqIO.parse(input_trinity_fasta_path, "fasta"))

	# Default dict that links the gene name to the largest variant size
	gene_representative_long_iso_dict = {}
	# Default dict that will hold the longest sequence record for each gene
	gene_current_max_iso_dict = {}

	print('Searching for longest isoforms')
	for i in range(len(trin_file)):
		seq_rec = trin_file[i]
		gene_id = '_'.join(seq_rec.id.split(' ')[0].split('_')[:4])
		gene_length = len(seq_rec.seq)
		if gene_id not in gene_current_max_iso_dict:  # Then no representative longest isoform so add this one
			gene_current_max_iso_dict[gene_id] = gene_length
			gene_representative_long_iso_dict[gene_id] = seq_rec
		else:  # Need to compare length of current isofrom to current longest isoform
			if gene_length > gene_current_max_iso_dict[gene_id]:  # If loner then this should be the new representative
				gene_current_max_iso_dict[gene_id] = gene_length
				gene_representative_long_iso_dict[gene_id] = seq_rec
			else:
				# the isoform of the gene was shorter than the current longest representative of the gene
				pass

	# Here the gene_representative_long_iso_dict contains the sequence records that we want to write out
	# Write this out manually so that we can avoid the automatic interleaving that seqIO does.
	print(f'Writing out new trinity assembly fasta to: {output_trinity_fasta_path}')
	with open(output_trinity_fasta_path, 'w') as f:
		for seq_rec in list(gene_representative_long_iso_dict.values()):
			f.write(f'>{seq_rec.description}\n')
			f.write(f'{seq_rec._seq}\n')
	# SeqIO.write(list(gene_representative_long_iso_dict.values()), output_trinity_fasta_path, "fasta")
	print('Done')

if __name__ == "__main__":
	remove_short_isoform_sequences()