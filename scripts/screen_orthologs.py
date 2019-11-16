import sys
import pandas as pd
from collections import defaultdict
import os
class screen_orthologs:
    def __init__(self):
        self.orth_df_in_single = pd.read_csv(filepath_or_buffer=sys.argv[1], header=0, sep='\t')
        self.orth_out_path = sys.argv[2]
        os.makedirs(os.path.dirname(self.orth_out_path), exist_ok=True)

    def screen_orfs(self):
        # This gave exactly the same result as for the orth_df_multi
        self.orth_df_in_single.set_index("group_id", drop=True, inplace=True)
        headers = list(self.orth_df_in_single)
        self.orth_df_in_single = self.orth_df_in_single[(self.orth_df_in_single[headers[0]] != '*') & (self.orth_df_in_single[headers[1]] != '*') & (self.orth_df_in_single[headers[2]] != '*') & (self.orth_df_in_single[headers[3]] != '*')]


        # Now we need to look to see that each of the transcipts is used only once
        # we can start by doing this for just the first strain
        species_headers = list(self.orth_df_in_single)
        for i in range(4):
            trans_used = defaultdict(list)
            too_many = []
            for ind, ser in self.orth_df_in_single.iterrows():
                trans_used[ser[species_headers[0]]].append(ser.name)
            # here we can see which of the transcripts have been used in more than one group
            for trans_name_key, orth_group_list in trans_used.items():
                if len(orth_group_list) > 1:
                    raise RuntimeError(f"A transciprt for {species_headers[0]} was used in {len(orth_group_list)} groups.")

        # It can seems that all of the transcipts were used in just one of the ortholog groupings.
        # now write out the dataframe
        self.orth_df_in_single.to_csv(path_or_buf=self.orth_out_path, sep='\t')


so = screen_orthologs()
so.screen_orfs()