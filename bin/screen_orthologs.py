#!/usr/bin/env python
import sys
import pandas as pd
from collections import defaultdict
import os
class screen_orthologs:
    def __init__(self):
        self.orth_df_in_single = pd.read_csv(filepath_or_buffer=sys.argv[1], header=0, sep='\t')
        self.orth_out_path = 'screened_orthologs.tsv'

    def screen_orfs(self):
        self.orth_df_in_single.set_index("group_id", drop=True, inplace=True)
        # The rows are ouput in an order such that the ortholog genes that have the most number of species in common
        # are put first. As such we can just go down the indices looking for an asterix. Once we find one
        # we know that we will need to drop all rows after that poipnt
        for i in range(len(self.orth_df_in_single.index.values.tolist())):
            if '*' in self.orth_df_in_single.iloc[i].values.tolist():
                ind_to_drop = i
                break
        self.orth_df_in_single = self.orth_df_in_single.iloc[:ind_to_drop]
        self.orth_df_in_single.to_csv(path_or_buf=self.orth_out_path, sep='\t')

so = screen_orthologs()
so.screen_orfs()