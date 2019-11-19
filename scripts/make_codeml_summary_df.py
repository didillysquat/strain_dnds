import sys
import os
import pandas as pd
import re
import pickle
import itertools
import numpy as np

class MakeCodemlSumDF:
    def __init__(self):
        self.species = sys.argv[1]
        self.df_pickle_out_path = sys.argv[2]
        self.df_csv_out_path = sys.argv[3]
        self.block_analysis_base_dir = os.path.join(
            '/home/humebc/projects/parky/breviolum_transcriptomes/codeml', self.species)
        self.list_of_block_dirs = self._populate_list_of_block_dirs()
        self.df_cols = self._make_df_cols()
        self.df = pd.DataFrame(columns=self.df_cols)
        self.count = 0
        self.dict_of_dnds_values = {}

    def _populate_list_of_block_dirs(self):
        list_of_dirs = list()
        for root, dirs, files in os.walk(self.block_analysis_base_dir):
            list_of_dirs = dirs
            break
        return list_of_dirs

    def _make_df_cols(self):
        cols = []
        if self.species == 'b_minutum':
            srr_list = ["SRR1793320", "SRR1793321", "SRR1793322", "SRR1793323"]
        else: #self.species == 'b_psygmophilum':
            srr_list = ["SRR1793324", "SRR1793325", "SRR1793326", "SRR1793327"]
        for srr_pair in itertools.combinations(srr_list, 2):
            cols.append(f'{srr_pair[0]}_{srr_pair[1]}_dn_ds')
            cols.append(f'{srr_pair[0]}_{srr_pair[1]}_dn')
            cols.append(f'{srr_pair[0]}_{srr_pair[1]}_ds')
        return cols


    def summarise_CODEML_pairwise_analyses(self):
        ''' Here we will go through each of the block analyses and take out the pairwise dn/ds analyses and put
        them into a single pandas dataframe. The ortholog will be the index (which we will evenutally sort) and the
        columns will be each of the pairwise comparisons.'''

        # the list that we will hold the info for a single row in
        for block_dir in self.list_of_block_dirs:
            sys.stdout.write('\n\nProcessing block {}\n\n'.format(block_dir))
            # read in the file for this block and populate the df with the dN/dS videos
            with open('{0}/{1}/{1}_codeml_results.out'.format(self.block_analysis_base_dir, block_dir), 'r') as f:
                out_file = [line.rstrip() for line in f]

            # go line by line through the output file picking up the dnds values and putting them into the df
            for i in range(len(out_file)):
                if 'Data set' in out_file[i]:
                    sys.stdout.write('\rProcessing {}'.format(out_file[i]))
                # see if the line in question contains dnds values
                if 'dN/dS=' in out_file[i]:
                    # check to see if we have populated a row worth of dnds values
                    # if so, populate the df and then start a new row list to collect values in
                    if self.count % 6 == 0 and self.count != 0:
                        # then we should have a set of dnds values that we can append to the df
                        holder_list = []
                        for col in self.df_cols:
                            holder_list.append(self.dict_of_dnds_values[col])
                        self.df.loc[int(ortholog_id)] = holder_list
                        dict_of_dnns_values = {}

                    # when we are here we are either adding one more value to an already exisiting set
                    # or we are starting a brand new set

                    # get the dn/ds value and the pair that we are working with
                    # num_reg_ex = re.compile('[0-9]+\.[0-9]+')
                    # num_reg_ex_d = re.compile('dN/dS\s*=\s*[0-9]+\.[0-9]+')
                    data_reg_ex = re.compile('dN/dS\s*=\s*([0-9]+\.[0-9]+)\s+dN\s*=\s*([0-9]+\.[0-9]+)\s+dS\s*=\s*([0-9]+\.[0-9]+)')
                    if 'nan' in out_file[i]:
                        dn_ds_value = np.nan
                        dn = np.nan
                        ds = np.nan
                    else:
                        data_match = data_reg_ex.findall(out_file[i])
                        dn_ds_value = float(data_match[0][0])
                        dn = float(data_match[0][1])
                        ds = float(data_match[0][2])
                    names_reg_ex = re.compile('\(([0-9]+)_(SRR[0-9]+)\).+\([0-9]+_(SRR[0-9]+)\)')
                    names_match = names_reg_ex.findall(out_file[i-4])

                    if ds > 60:
                        apples = 'asdf'
                    ortholog_id = int(names_match[0][0])
                    spp_one = names_match[0][1]
                    spp_two = names_match[0][2]

                    if '{}_{}_dn_ds'.format(spp_one, spp_two) in self.df_cols:
                        self.dict_of_dnds_values['{}_{}_dn_ds'.format(spp_one, spp_two)] = dn_ds_value
                        self.dict_of_dnds_values['{}_{}_dn'.format(spp_one, spp_two)] = dn
                        self.dict_of_dnds_values['{}_{}_ds'.format(spp_one, spp_two)] = ds
                    else:
                        self.dict_of_dnds_values['{}_{}_dn_ds'.format(spp_two, spp_one)] = dn_ds_value
                        self.dict_of_dnds_values['{}_{}_dn'.format(spp_two, spp_one)] = dn
                        self.dict_of_dnds_values['{}_{}_ds'.format(spp_two, spp_one)] = ds

                    self.count += 1

        # # finally add the last set of dn/ds data to the df
        # holder_list = []
        # for col in self.df_cols:
        #     holder_list.append(dict_of_dnns_values[col])
        # self.df.loc[int(ortholog_id)] = holder_list


        pickle.dump(self.df, open(self.df_pickle_out_path, 'wb'))
        self.df.to_csv(self.df_csv_out_path)

mcsdf = MakeCodemlSumDF()
mcsdf.summarise_CODEML_pairwise_analyses()