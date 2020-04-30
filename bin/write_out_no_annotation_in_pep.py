#!/usr/bin/env python3
# This script will take in the output of the panzer annotation files from using the
# Panzer2 server.
# We will aim to add to the already made go df for each transciptome
# We will also aim to output a final .pep file that once again contains the no go annotation
# files. We will also pickle out another dict.

import sys
import os
from collections import defaultdict
import pandas as pd
import pickle
import numpy as np

class WriteOutNoAnnoInFasta:
    def __init__(self):
        # read in the old pep file
        self.pep_input_directory = sys.argv[2]
        self.original_pep_file_name = sys.argv[3].split('/')[-1]
        try:
            with open(os.path.join(self.pep_input_directory, self.original_pep_file_name), 'r') as f:
                pep_file_list = [line.rstrip() for line in f]
                self.pep_file_dict = {pep_file_list[i].split()[0][1:]:pep_file_list[i+1] for i in range(0, len(pep_file_list), 2)}
        except FileNotFoundError as e:
            print(e)
            raise RuntimeError("Couldn't locate original .pep file")
        
        # Read in the current go_df_dict
        self.base_name = self.original_pep_file_name.replace('.mmseqs.panzer.in.pep', '')
        self.cache_dir_base = sys.argv[4]
        self.go_df_dict_path = os.path.join(self.cache_dir_base, f'{self.base_name}_go_df_dict.p')
        self.go_df_dict = pickle.load(open(self.go_df_dict_path, 'rb'))
        self.go_df_dict_w_panzer_path = f'{self.base_name}_go_df_w_panzer_dict.p'
        
        # Read in the panzer.go.out
        
        self.panzer_go_out_path = os.path.join(sys.argv[1], f'{self.base_name}.panzer.go.out')
        with open(self.panzer_go_out_path, 'r') as f:
            panzer_list = [_.rstrip() for _ in f]
        
        # Let's make a quick dict to see how many queries we are dealing with
        query_to_go_dict = defaultdict(set)
        for line in panzer_list[1:]:
            query = line.split('\t')[0]
            go_term = line.split('\t')[2]
            query_to_go_dict[query].add(f"GO:{go_term}")
        
        # Now for each of the queries have a quick look to see if it is already in
        # the go df dict. If not then put it in using nan for the evalue
        count_in = 0
        for k, v in query_to_go_dict.items():
            if k not in self.go_df_dict:
                self.go_df_dict[k] = ['panzer', np.nan, np.nan, np.nan, ';'.join(list(v))]
            else:
                count_in += 1
        print(f"{count_in} were already in go_df_dict")
        self.new_pep_file_name = self.original_pep_file_name.replace('.mmseqs.panzer.in.pep', '.mmseqs.no_annotation.in.pep')
        seqs_for_pep = [k for k in self.pep_file_dict if k not in self.go_df_dict]
        self.new_pep_list = []
        for query_name in seqs_for_pep:
            self.new_pep_list.extend([f'>{query_name}', self.pep_file_dict[query_name]])
        with open(self.new_pep_file_name, 'w') as f:
            f.write('\n'.join(self.new_pep_list))
        
        # pickle out the new dict too
        pickle.dump(self.go_df_dict, open(self.go_df_dict_w_panzer_path, 'wb'))

        # create the df and write out
        df = pd.DataFrame.from_dict(self.go_df_dict, orient='index', columns=['db', 'match_accession', 'e_value', 'match_ranking', 'go_terms'])
        df.index.name = 'query'
        df.to_csv(f'{self.base_name}_go_df_w_panzer.csv', index=True)
        

if __name__ == "__main__":
    wotif = WriteOutNoAnnoInFasta()