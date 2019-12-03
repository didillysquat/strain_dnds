import sys
import pandas as pd
import itertools
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import os

class ProcessData:
    def __init__(self):
        self.min_df = pd.read_csv('/home/humebc/projects/parky/breviolum_transcriptomes/codeml/b_minutum/codeml_summary_df.csv', index_col=0)
        self.psy_df = pd.read_csv(
            '/home/humebc/projects/parky/breviolum_transcriptomes/codeml/b_psygmophilum/codeml_summary_df.csv', index_col=0)
        self.species_list = ['b_minutum', 'b_psygmophilum']
        self.sra_dict = {"b_minutum":["SRR1793320", "SRR1793321", "SRR1793322", "SRR1793323"], "b_psygmophilum": ["SRR1793324", "SRR1793325", "SRR1793326", "SRR1793327"]}
        self.fig_dir = '/home/humebc/projects/parky/breviolum_transcriptomes/figures'

    def explore(self):
        foo = 'bar'
        # first lets try to plot up how many psgs we have going on

        fig, ax_arr = plt.subplots(nrows=2)
        for species in self.species_list:
            if species == 'b_minutum':
                df = self.min_df
                ax = ax_arr[0]
            else:
                df = self.psy_df
                ax = ax_arr[1]
            sra = self.sra_dict[species]
            plot_dict = defaultdict(list)
            for cut_off in np.arange(0, 0.1, 0.001):
                for sp_1, sp_2 in itertools.combinations(sra, 2):
                    if f'{sp_1}_{sp_2}_dn_ds' in list(df):
                        base_col_name = f'{sp_1}_{sp_2}'
                    else:
                        base_col_name = f'{sp_2}_{sp_1}'

                    # for each cutoff we want to get a pari of psgs and cutoff points for each of the species
                    psg_df = df[(df[f'{base_col_name}_dn_ds'] > 1) & (df[f'{base_col_name}_ds'] > cut_off)]
                    num_psg = len(psg_df.index)
                    plot_dict[base_col_name].append(num_psg)

            #here we have the data to be plotted
            plot_df = pd.DataFrame.from_dict(plot_dict, orient='index', columns=list(np.arange(0, 0.1, 0.001)))
            plot_df = plot_df.T
            plot_df['x'] = plot_df.index.values.tolist()

            ax.set_xlim(0,0.01)
            ax.set_ylim(0,300)
            for comparison in list(plot_df)[:-1]:
                ax.plot('x', comparison, data=plot_df, marker='', color='black', linewidth=2)
            foo = 'bar'
        plt.savefig(os.path.join(self.fig_dir, 'dn_ds_with_ds_cutoff_both_species.pdf'))


pdata = ProcessData()
pdata.explore()