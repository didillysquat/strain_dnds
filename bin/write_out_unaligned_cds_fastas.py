#!/usr/bin/env python
"""
Use the screened_orthologs.tsv to pull out the ortholog genese from the .cds
set of files. For each gene, write out a fasta that contains the gene representative
for each of the strains.
The column titles of the screened ortholog .tsv are the .pep file
names (e.g. SRR1793320_longest_iso_orfs.single_orf.pep)
The cds files that we will be getting the sequences from are all prefixed
with the SRRXXX (e.g. SRR1793320_longest_iso_orfs.cds). As such we can use
this SRRXXX in the python script to map between the table and the .cds files"""
import sys
import os
import pandas as pd
from Bio import SeqIO

class WriteOutFastas:
    def __init__(self):
        self.path_to_ortholog_table = os.path.abspath(sys.argv[1])
        self.orth_df = pd.read_csv(self.path_to_ortholog_table, header=0, sep='\t', index_col=0)
        # rename the headers so that they are only the strain IDs
        # They are currently e.g. SRR1793320_longest_iso_orfs.single_orf.pep
        # This will make them SRR1793320
        self.strain_names = [head.replace('_longest_iso_orfs.single_orf.pep', '') for head in list(self.orth_df)]
        self.orth_df.columns = self.strain_names
        # Make a dict so that we can quickly grab the sequence (nucleotide not protein) of a transcipt in question
        # it will be a nested dict so that first key is strain and second it transcript name
        self.transcript_name_to_seq_dict = {}
        self._init_transcript_name_to_seq_dict()
        self.non_three_orth_groups = []

    def _init_transcript_name_to_seq_dict(self):
        print("Creating transcript name to cds sequence mapping")
        for strain_name in self.strain_names:
            print(strain_name)
            transcript_name_to_seq_dict_temp = {}
            # These files should all be in the working directory when we are running this on nextflow
            cds_file_name = f'{strain_name}_longest_iso_orfs.cds'
            strain_cds_fasta_as_file = list(SeqIO.parse(cds_file_name, "fasta"))
            for seq_record in strain_cds_fasta_as_file:
                transcript_name_to_seq_dict_temp[seq_record.id] = str(seq_record.seq)
            self.transcript_name_to_seq_dict[strain_name] = transcript_name_to_seq_dict_temp

    def write_out_fasta(self):
        # We will work through the orth_df and write
        # 090320 I am modifying this so that we write the files out directly to the directory
        # rather than having a new subdirectory for each gene
        print("writing out local unaligned fastas")
        for ind, ser in self.orth_df.iterrows():
            sys.stdout.write(f"\r{ind}")
            orth_group_unaligned_fasta_path = f'{ind}_unaligned_cds.fasta'
            fasta = []
            non_three = False
            for strain_name in self.strain_names:
                # check to make sure that the seq is divisble by three
                seq_str = self.transcript_name_to_seq_dict[strain_name][ser[strain_name]]
                if (len(seq_str)%3) != 0:
                    non_three = True
                fasta.extend([f'>{ind}_{strain_name}', seq_str])
            if non_three:
                self.non_three_orth_groups.append(str(ind))
                continue
            with open(orth_group_unaligned_fasta_path, 'w') as f:
                for line in fasta:
                    f.write(f'{line}\n')

wof = WriteOutFastas()
wof.write_out_fasta()