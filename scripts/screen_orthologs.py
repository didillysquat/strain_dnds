import sys
import pandas as pd
from collections import defaultdict
class screen_orthologs:
    def __init__(self):
        self.orth_df_in_multi = pd.read_csv(filepath_or_buffer=sys.argv[1], header=0, sep='\t')
        self.orth_df_in_single = pd.read_csv(filepath_or_buffer=sys.argv[2], header=0, sep='\t')
        self.orth_out_path = sys.argv[3]

    def screen_orfs(self):
        # only keep the orthologs that were found in all 8 of the strains
        self.orth_df_in_multi = self.orth_df_in_multi[self.orth_df_in_multi['sp_in_grp'] == 8]
        self.orth_df_in_multi = self.orth_df_in_multi[self.orth_df_in_multi['group_size'] == 8]
        # only keep the rows that do not contain *
        self.orth_df_in_single.set_index(keys='group_id', drop=True, inplace=True)
        print(f'The index is {len(self.orth_df_in_single.index)} long')
        for col in list(self.orth_df_in_single):
            self.orth_df_in_single = self.orth_df_in_single[self.orth_df_in_single[col] != '*']
            print(f'The index is {len(self.orth_df_in_single.index)} long')

        foo = 'bar'

        # Now we need to look to see that each of the transcipts is used only once
        # we can start by doing this for just the first strain
        trans_used = defaultdict(list)
        for ind, ser in self.orth_df_in_single.iterrows():
            for r_ind in ser.index.values.tolist():
                trans_used[ser[r_ind]].append(ser.name)

        foo = "bar"

        # First lets break it down to only those transcripts that have exact matches in each of the gen

        # Dictionary to keep track of which of the transcripts have been used


so = screen_orthologs()
so.screen_orfs()