import os
with open('/home/humebc/projects/parky/breviolum_transcriptomes/nf_master_fasta_and_q_file/master_fasta_for_tree.fasta', 'r') as f:
    fasta_file = [line.rstrip() for line in f]

seqs = [fasta_file[i+1] for i in range(0,len(fasta_file), 2)]
num_tax_range = range(len(seqs))
yes_count = 0
for i in range(len(seqs[0])):
    aa_set = set([seqs[j][i] for j in num_tax_range])
    if (len(aa_set) > 1) and ('-' not in aa_set):
        yes_count += 1
print(f'{yes_count} informative sites out of {len(seqs[0])}')
    