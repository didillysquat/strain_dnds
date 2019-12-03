import sys
import os
import pandas as pd
from Bio import SeqIO

class WriteOutFastas:
    def __init__(self):
        self.species = sys.argv[1]
        self.path_to_ortholog_table = os.path.abspath(sys.argv[2])
        self.orf_prediction_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(sys.argv[3])), '..'))
        self.orth_df = pd.read_csv(self.path_to_ortholog_table, header=0, sep='\t', index_col=0)
        # rename the headers so that they are only the strain IDs
        self.strain_names = [head.split('_')[0] for head in list(self.orth_df)]
        self.orth_df.columns = self.strain_names
        self.output_path = os.path.abspath(os.path.join(os.path.dirname(self.path_to_ortholog_table), '..', '..', "local_alignments", self.species))
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
            cds_path = os.path.join(self.orf_prediction_dir, strain_name, 'longest_iso_orfs.cds')
            strain_cds_fasta_as_file = list(SeqIO.parse(cds_path, "fasta"))
            for seq_record in strain_cds_fasta_as_file:
                transcript_name_to_seq_dict_temp[seq_record.id] = str(seq_record.seq)
            self.transcript_name_to_seq_dict[strain_name] = transcript_name_to_seq_dict_temp

    def write_out_fasta(self):
        # We will work through the orth_df and write
        print("writing out local unaligned local fastas")

        for ind, ser in self.orth_df.iterrows():
            sys.stdout.write(f"\r{ind}")
            orth_group_output_dir = os.path.join(self.output_path, str(ind))
            orth_group_unaligned_fasta_path = os.path.join(orth_group_output_dir, f'{ind}_unaligned_cds.fasta')
            if ind == 12175:
                foo = 'bar'
            fasta = []
            non_three = False
            for strain_name in self.strain_names:
                # check to make sure that the seq is divisble by three
                seq_str = self.transcript_name_to_seq_dict[strain_name][ser[strain_name]]
                if len(seq_str) == 0:
                    foo = 'bar'
                if (len(seq_str)%3) != 0:
                    non_three = True
                fasta.extend([f'>{ind}_{strain_name}', seq_str])
            if non_three:
                self.non_three_orth_groups.append(str(ind))
                continue
            os.makedirs(orth_group_output_dir, exist_ok=True)
            with open(orth_group_unaligned_fasta_path, 'w') as f:
                for line in fasta:
                    f.write(f'{line}\n')
        
        self._write_out_summary_file()

    def _write_out_summary_file(self):
        # We will write out a summary file so that we have an output for the snakemake to grab hold of
        print("Writing out summary file")
        with open(os.path.join(self.output_path, "unaligned_fastas_summary.readme"), 'w') as f:
            f.write(f"There were {len(self.non_three_orth_groups)} groups containing at least one "
                    f" sequence that was not divisible by three\n")

wof = WriteOutFastas()
wof.write_out_fasta()