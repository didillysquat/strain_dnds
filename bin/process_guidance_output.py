#!/usr/bin/env python
"""
Perform a mafft alignment from inside the guidance software 
that will produce a cropped and aligned cds and aa alignment
We will incorporate dropping columns of low support and 
removing parts of the alignment at the 
beginning and end for which any one of the orthologs has nucleotides missing.
This script will work on a single unaligned cds fasta.
It will output the two output the aligned cds and aa as:
xxx_cropped_aligned_aa.fasta  xxx_cropped_aligned_cds.fasta
"""
import os
import sys
from multiprocessing import Queue, Process
import ntpath
import subprocess
import pandas as pd
import shutil
import ntpath

class ProcessGuidanceOutput:
    def __init__(self):
        self.seq_names = None
        self.orth_group_id = ntpath.basename(sys.argv[1]).split('.')[0]
        print(f'orth_group_id is {self.orth_group_id}')
        self.aa_cols_score_file_path = sys.argv[1]
        self.aa_alignment_file_path = sys.argv[2]
        # This aa_alignment has the original names in it. 
        # we use it to get the names in the _get_seq_names_from_fasta method below
        self.cds_alignment_file_path = sys.argv[3]
        self._get_seq_names_from_fasta()
        # Dataframes that will hold the aa and cds alignments that will be used for cropping and col dropping
        self.aa_df = None
        self.cds_df = None
        self.aa_cols_to_remove = []
        self.cds_cols_to_remove = []
        # output paths for the aligned and cropped aa and cds alignments
        self.aa_aligned_and_cropped_path = f'{self.orth_group_id}_cropped_aligned_aa.fasta'
        self.cds_aligned_and_cropped_path = f'{self.orth_group_id}_cropped_aligned_cds.fasta'

    def _get_seq_names_from_fasta(self):
        with open(self.cds_alignment_file_path, 'r') as f:
            fasta_file = self.convert_interleaved_to_sequencial_fasta_two([line.rstrip() for line in f])
        self.seq_names = [fasta_file[i][1:] for i in range(0,len(fasta_file),2)]

    def drop_low_qual_cols_and_crop_and_write(self):
        """
        At this point we have performed the guidance analysis.
        We now read in the outputs that tell us the scores of the aa and cds alignment columns
        we then use these to drop columns that have a low score.
        We also then crop the alignment in the begining and the end so that
        any columns for which one of th sequences did not have a sequence are dropped.
        """

        self._get_aa_cols_to_remove()

        # here we have a 0 based list of the columns that need to be dropped from the aa alignment
        # convert this to a 0 based index of the columns that need to be dropped from the cds alignment
        self._get_cds_cols_to_remove()

        # here we have the indices that need to be dropped for the cds and aa alignments
        # now read in the alignments as 2D lists, convert to pandas dataframe and perform the columns drops

        # aa first
        self._make_aa_df()

        # now drop the columns
        self._drop_cols_from_aa_df()

        # cds second
        self._make_cds_df()

        # now drop the columns
        self._drop_cols_from_cds_df()

        # here we have the cds and aa dfs that we can now do the cropping with and then finally write back out as
        # fasta files
        # go from either end deleting any columns that have a gap

        # aa first
        self.aa_df = self._crop_fasta_df(aligned_fasta_as_pandas_df_to_crop=self.aa_df)

        # cds second
        self.cds_df = self._crop_fasta_df(aligned_fasta_as_pandas_df_to_crop=self.cds_df)

        # now we just need to write out dfs
        self._write_out_aa_df()

        self._write_out_cds_df()

        # here we have the cds and the aa alignments cropped and written out
        # we can then use these as input into CODEML and the BUSTED programs

    def _write_out_cds_df(self):
        cds_fasta_list = self.pandas_df_to_fasta(pd_df=self.cds_df, seq_names=self.seq_names)
        with open(self.cds_aligned_and_cropped_path, 'w') as f:
            for line in cds_fasta_list:
                f.write(f'{line}\n')

    def _write_out_aa_df(self):
        aa_fasta_list = self.pandas_df_to_fasta(pd_df=self.aa_df, seq_names=self.seq_names)
        with open(self.aa_aligned_and_cropped_path, 'w') as f:
            for line in aa_fasta_list:
                f.write(f'{line}\n')

    def _drop_cols_from_cds_df(self):
        columns_to_keep_cds = [col for col in list(self.cds_df) if col not in self.cds_cols_to_remove]
        self.cds_df = self.cds_df[columns_to_keep_cds]

    def _make_cds_df(self):
        with open(self.cds_alignment_file_path, 'r') as f:
            cds_alignment_file_list = self.convert_interleaved_to_sequencial_fasta_two([line.rstrip() for line in f])
        self.cds_df = pd.DataFrame([list(line) for line in cds_alignment_file_list if not line.startswith('>')])

    def _drop_cols_from_aa_df(self):
        columns_to_keep_aa = [col for col in list(self.aa_df) if col not in self.aa_cols_to_remove]
        self.aa_df = self.aa_df[columns_to_keep_aa]

    def _make_aa_df(self):
        with open(self.aa_alignment_file_path, 'r') as f:
            aa_alignment_file_list = self.convert_interleaved_to_sequencial_fasta_two([line.rstrip() for line in f])
        self.aa_df = pd.DataFrame([list(line) for line in aa_alignment_file_list if not line.startswith('>')])

    def _get_cds_cols_to_remove(self):
        for col_index in self.aa_cols_to_remove:
            self.cds_cols_to_remove.extend([col_index * 3 + n for n in range(3)])

    def _get_aa_cols_to_remove(self):
        with open(self.aa_cols_score_file_path, 'r') as f:
            aa_col_score_file_list = [line.rstrip() for line in f][1:]
        for i in range(int(len(aa_col_score_file_list) / len(self.seq_names))):
            # the range that represents i and the next three i s in the iteration without increasing i
            file_line_indices = [i * len(self.seq_names) + k for k in range(len(self.seq_names))]
            scores_set = [
                float(aa_col_score_file_list[j].split()[2]) if '-nan' not in aa_col_score_file_list[j] else '-nan' for j
                in file_line_indices]

            # now examine what is in the score sets
            drop = False
            nan_count = 0
            for score in scores_set:
                if score != '-nan':
                    if score < 0.6:
                        drop = True
                        break
                elif score == '-nan':
                    nan_count += 1
                else:
                    continue
            # when we come out ned to check if the nan_score == 4
            if nan_count == len(self.seq_names):
                drop = True

            # if drop = True then this is a column to drop
            if drop:
                self.aa_cols_to_remove.append(i)

    def convert_interleaved_to_sequencial_fasta_two(self, fasta_in):
        fasta_out = []

        for i in range(len(fasta_in)):

            if fasta_in[i].startswith('>'):
                if fasta_out:
                    # if the fasta is not empty then this is not the first
                    fasta_out.append(temp_seq_str)
                # else then this is the first sequence and there is no need to add the seq.
                temp_seq_str = ''
                fasta_out.append(fasta_in[i])
            else:
                temp_seq_str = temp_seq_str + fasta_in[i]
        # finally we need to add in the last sequence
        fasta_out.append(temp_seq_str)
        return fasta_out

    def _crop_fasta_df(self, aligned_fasta_as_pandas_df_to_crop):
        columns_to_drop = []
        for i in list(aligned_fasta_as_pandas_df_to_crop):
            # if there is a gap in the column at the beginning
            if '-' in list(aligned_fasta_as_pandas_df_to_crop[i]) or '*' in list(aligned_fasta_as_pandas_df_to_crop[i]):
                columns_to_drop.append(i)
            else:
                break
        for i in reversed(list(aligned_fasta_as_pandas_df_to_crop)):
            # if there is a gap in the column at the end
            if '-' in list(aligned_fasta_as_pandas_df_to_crop[i]) or '*' in list(aligned_fasta_as_pandas_df_to_crop[i]):
                columns_to_drop.append(i)
            else:
                break

        # get a list that is the columns indices that we want to keep
        col_to_keep = [col_index for col_index in list(aligned_fasta_as_pandas_df_to_crop) if
                       col_index not in columns_to_drop]
        # drop the gap columns
        return aligned_fasta_as_pandas_df_to_crop[col_to_keep]

    def pandas_df_to_fasta(self, pd_df, seq_names):
        temp_fasta = []
        for i in range(len(seq_names)):
            temp_fasta.extend(['>{}'.format(seq_names[i]), ''.join(list(pd_df.loc[i]))])
        return temp_fasta

pgo = ProcessGuidanceOutput()
pgo.drop_low_qual_cols_and_crop_and_write()