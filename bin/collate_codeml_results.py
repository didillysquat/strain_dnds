#!/usr/bin/env python3
import sys
import os
import pandas as pd
import re
import itertools
import numpy as np

class MakeCodemlSumDF:
    def __init__(self):
        self.strain_list = sys.argv[1].split(',')
        self.codeml_out_path_list = self._populate_codeml_out_path_list()
        self.df_cols = self._make_df_cols()
        self.df = pd.DataFrame(columns=self.df_cols)
        self.count = 0
        self.dict_of_dnds_values = {}
        self.num_combinations = len(list(itertools.combinations(self.strain_list, 2)))
        self.data_reg_ex = re.compile('dN/dS\s*=\s*([0-9]+\.[0-9]+)\s+dN\s*=\s*([0-9]+\.[0-9]+)\s+dS\s*=\s*([0-9]+\.[0-9]+)')
        self.names_reg_ex = re.compile('\(([0-9]+)_(\S+)\).+\([0-9]+_(\S+)\)')

    def _populate_codeml_out_path_list(self):
        files = list(os.walk('.'))[0][2]
        return [file_name for file_name in files if '.out' in file_name]
        
    def _make_df_cols(self):
        cols = []
        for srr_pair in itertools.combinations(self.strain_list, 2):
            cols.append(f'{srr_pair[0]}_{srr_pair[1]}_dn_ds')
            cols.append(f'{srr_pair[0]}_{srr_pair[1]}_dn')
            cols.append(f'{srr_pair[0]}_{srr_pair[1]}_ds')
        return cols

    def summarise_CODEML_pairwise_analyses(self):
        ''' Here we will go through each of codeml outputs and take out the pairwise dn/ds analyses and put
        them into a single pandas dataframe. The ortholog will be the index (which we will evenutally sort) and the
        columns will be each of the pairwise comparisons.'''

        # the list that we will hold the info for a single row in
        for codeml_out_file in self.codeml_out_path_list:
            self.dict_of_dnds_values = {}
            sys.stdout.write(f'\r{codeml_out_file}')
            # read in the file for this block and populate the df with the dN/dS values
            with open(codeml_out_file, 'r') as f:
                out_file = [line.rstrip() for line in f]

            # go line by line through the output file picking up the dnds values and putting them into the df
            for i in range(len(out_file)):
                if 'Data set' in out_file[i]:
                    sys.stdout.write('\rProcessing {}'.format(out_file[i]))
                # see if the line in question contains dnds values
                if 'dN/dS=' in out_file[i]:
                    # get the dn/ds value and the pair that we are working with
                    if 'nan' in out_file[i]:
                        dn_ds_value = np.nan
                        dn = np.nan
                        ds = np.nan
                    else:
                        data_match = self.data_reg_ex.findall(out_file[i])
                        dn_ds_value = float(data_match[0][0])
                        dn = float(data_match[0][1])
                        ds = float(data_match[0][2])

                    names_match = self.names_reg_ex.findall(out_file[i-4])
                    ortholog_id = int(names_match[0][0])  
                    spp_one = names_match[0][1]
                    spp_two = names_match[0][2]

                    if '{}_{}_dn_ds'.format(spp_one, spp_two) in self.df_cols:
                        self.dict_of_dnds_values['{}_{}_dn_ds'.format(spp_one, spp_two)] = dn_ds_value
                        self.dict_of_dnds_values['{}_{}_dn'.format(spp_one, spp_two)] = dn
                        self.dict_of_dnds_values['{}_{}_ds'.format(spp_one, spp_two)] = ds
                        self.count += 1
                    elif '{}_{}_dn_ds'.format(spp_two, spp_one) in self.df_cols:
                        self.dict_of_dnds_values['{}_{}_dn_ds'.format(spp_two, spp_one)] = dn_ds_value
                        self.dict_of_dnds_values['{}_{}_dn'.format(spp_two, spp_one)] = dn
                        self.dict_of_dnds_values['{}_{}_ds'.format(spp_two, spp_one)] = ds
                        self.count += 1
                    else:
                        raise RuntimeError(f'neither {spp_one}_{spp_two}_dn_ds nor {spp_two}_{spp_one}_dn_ds were found in the df columns')

                    # check to see if we have populated a row worth of dnds values
                    # if so, populate the df and then start a new row list to collect values in
                    if self.count % self.num_combinations == 0 and self.count != 0:
                        # then we should have a set of dnds values that we can append to the df
                        holder_list = []
                        for col in self.df_cols:
                            holder_list.append(self.dict_of_dnds_values[col])
                        self.df.loc[int(ortholog_id)] = holder_list
                        # here we now need to move on to the next codeml_out_file
                        self.count = 0
                        break
            
        self.df.to_csv("codeml_results_df.csv")

mcsdf = MakeCodemlSumDF()
mcsdf.summarise_CODEML_pairwise_analyses()