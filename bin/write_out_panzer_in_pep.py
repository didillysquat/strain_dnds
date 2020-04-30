#!/usr/bin/env python3
# This script will take in the output mmseq files from running the nr search
# and look to see which of the matches allow us to associate a GO.
# If we are not able to associate a GO then we will write out the query sequence
# into a new pep that will be run against the panzer server.
# For each query that we do get a hit for, we will add it to the go_df_dict
# that was created and pickled out as part of the processing of the sprot out
# results.
# To do the first mapping we have created a mapping dict that is pickled out
# at /home/humebc/projects/parky/strain_dnds/nf_gene2accession_mapping/gene2accession.p
# We can read this in. If a match is not in the keys, then there was no translation to
# a gene id.
# we will then need to read in the gene2go and make a dict from this.
# We are expecting the return from this to be very low.
# We should also convert the go_df_dict into a final dataframe that can be output
# for john.

import sys
import os
from collections import defaultdict
import pandas as pd
import pickle

class WriteOutPanzerInFasta:
    def __init__(self):
        # A set of the seq names that returned a result
        # And the actual output match file as a list
        self.mmseqs_out_name_set, self.mmseqs_match_list = self._make_mmseqs_out_name_set()
        # List of the full paths in which to search for the pep file of interest
        self.pep_input_directory = sys.argv[2]
        self.original_pep_file_name = sys.argv[3].split('/')[-1]
        self.cache_dir_base = sys.argv[4]
        
        self.base_name = self.original_pep_file_name.replace('.mmseqs.nr.in.pep', '')
        self.pep_file_dict = None
        self.df_columns = []
        self.go_count = 0
        try:
            with open(os.path.join(self.pep_input_directory, self.original_pep_file_name), 'r') as f:
                pep_file_list = [line.rstrip() for line in f]
                self.pep_file_dict = {pep_file_list[i].split()[0][1:]:pep_file_list[i+1] for i in range(0, len(pep_file_list), 2)}
        except FileNotFoundError as e:
            print(e)
            raise RuntimeError("Couldn't locate original .pep file")
        
        # the df that will hold the GO annotations that we were able to find.
        # We will make this as a dict to start with and then create the df in the final
        # script.
        self.go_df_dict_path = os.path.join(self.cache_dir_base, f'{self.base_name}_go_df_dict.p')
        self.go_df_dict = pickle.load(open(self.go_df_dict_path, 'rb'))
        # Read through the original .pep file and see which 
        # sequences didn't get a blast results
        # put these into the new pep
        self.new_pep_list = self._add_no_match_seqs_to_pep()
        
        # Make the look up dicts for gene2accession
        self.gene2accession_dict_pickle_path = sys.argv[5]
        # make sure it is no longer a default dict
        self.g2a_dict = dict(pickle.load(open(self.gene2accession_dict_pickle_path, 'rb')))
        self.gene2go_path = sys.argv[6]
        self.g2go_dict = self._make_g2go_dict()

        # now we need to go through the mmseqs file and see which matches can be associated to a GO
        self._log_go_associations()

        self.new_pep_file_name = self.original_pep_file_name.replace('.mmseqs.nr.in.pep', '.mmseqs.panzer.in.pep')
        self._write_out_new_pep()
        # Pickle out the go_df_dict
        pickle.dump(self.go_df_dict, open(self.go_df_dict_path, 'wb'))
        # But also write out the dict as a df finally
        df = pd.DataFrame.from_dict(self.go_df_dict, orient='index', columns=['db', 'match_accession', 'e_value', 'match_ranking', 'go_terms'])
        df.index.name = 'query'
        df.to_csv(f'{self.base_name}_go_df_dict.csv', index=True)

        
    def _make_g2go_dict(self):
        g2go_dict = defaultdict(set)
        with open(self.gene2go_path, 'r') as f:
            for line in f.readlines()[1:]:
                g2go_dict[line.split('\t')[1]].add(line.split('\t')[2])
        # We can slim this down to only include the geneids that are in the g2a dict
        gene_ids = set()
        for v in self.g2a_dict.values():
            gene_ids.update(v)
        g2go_dict = {k: v for k, v in g2go_dict.items() if k in gene_ids}
        # convert back from a default dict
        return g2go_dict
    
    def _log_go_associations(self):
        """
        This will do a lot of the work.
        We will work our way through the self.mmseqs_match_list
        We will process one match per query at a time.
        In this way, we will minimise the number of remote requests to the uniprot server
        we have to do.
        I'm worried that if we hit it with too big a request that we will cuase it some issues
        We will keep track of which result number we are on as we will want to log this.
        """
        match_number = 1
        # this will be a sub set of the results from the self.mmseqs_match_list
        # one for each query that hasn't already been associated with a GO
        # sub_match_list = []
        # sub_match_query_set = set()
        # sub_match_list_for_next_iteration = []
        while True:
            sub_match_list = []
            sub_match_query_set = set()
            sub_match_list_for_next_iteration = []
            for line in self.mmseqs_match_list:
                query = line.split('\t')[0]
                if query in self.go_df_dict:
                    # Then we already have an annotation for this query
                    # Simply skip this line
                    continue
                if query not in sub_match_query_set:
                    # Then we haven't had a GO annotation for this yet
                    # and we haven't already added one match for this query to the sub_match_list
                    # Add this match, add query to set
                    sub_match_list.append(line)
                    sub_match_query_set.add(query)
                elif query in sub_match_query_set:
                    # Then we already have a match to process for this query
                    # so we should add all other lines into the sub_match_list_for_next_iteration
                    sub_match_list_for_next_iteration.append(line)
            # At this point we have a list of query matches (one per query)
            # That we need to see if there are annotations for.
            # First we need to get a dictionary that maps the match accesssion to a list
            # of GO terms.
            print(f'\n{len(sub_match_query_set)} queries left to find GO matches for')
            print(f'{len(sub_match_list_for_next_iteration)} match lines left to process')
            
            self._associate_go(sub_match_list, match_number)
            
            if not sub_match_list_for_next_iteration:
                break
            else:
                self.mmseqs_match_list = sub_match_list_for_next_iteration
            match_number += 1
        
        # At this point we have found all of the go associations that we can from
        # this set of mmseqs search results
        # We need to create and pickle out the annotation dict
        # We need to add those queries that didn't return a valid annotation to the pep
        # for writing out for trembl.
        name_to_add_to_pep = [_ for _ in self.mmseqs_out_name_set if _ not in self.go_df_dict]
        for query_name in name_to_add_to_pep:
            self.new_pep_list.extend([f'>{query_name}', self.pep_file_dict[query_name]])
    
    def _associate_go(self, sub_match_list, match_number):
        # Now we have the mappings that we need to see if there are go associations
            for line in sub_match_list:
                match = line.split('\t')[1]
                # The match will only be in the dict if there is an annotation for it
                # log the association
                try:
                    # will return a set
                    gene_id = self.g2a_dict[match]
                    go_set = set()
                    for g_id in gene_id:
                        try:
                            go_set.update(self.g2go_dict[g_id])
                        except KeyError:
                            pass
                    if go_set:
                        self.go_count += 1
                        print(f'Now found {self.go_count}')
                        self.go_df_dict[line.split('\t')[0]] = ['nr', match, line.split('\t')[2], match_number, ';'.join(list(go_set))]
                except KeyError:
                    pass

    def _make_mmseqs_out_name_set(self):
        with open(sys.argv[1], 'r') as f:
            mmseqs_out_list = [line.rstrip() for line in f]
        
        mmseqs_out_name_set = set()
        for line in mmseqs_out_list:
            mmseqs_out_name_set.add(line.split('\t')[0])
        return mmseqs_out_name_set, mmseqs_out_list

    def _write_out_new_pep(self):
        with open(self.new_pep_file_name, 'w') as f:
            for line in self.new_pep_list:
                f.write(f'{line}\n')
        
    def _add_no_match_seqs_to_pep(self):
        """
        Work through the original pep file to see which seqs returned a blast result
        For those that didn't, add them to the new pep so that they can be blasted against
        the trembl db.
        """
        new_pep_list = []
        for k, v in self.pep_file_dict.items():
            if k not in self.mmseqs_out_name_set:
                # Then we had no blast result for this seq and it needs to be added to the
                # new pep file
                new_pep_list.extend([f'>{k}', v])
            else:
                pass
        return new_pep_list

# class CheekyFix:
#     def __init__(self):
#         self.list_of_dict_paths = self._make_list_of_dict_paths()
#         for path in self.list_of_dict_paths:
#             base_name = path.split('/')[-1].replace('_go_df_dict.p', '')
#             dict_to_fix = pickle.load(open(path, 'rb'))
#             path_to_trembl_in_pep = os.path.join("/home/humebc/projects/parky/strain_dnds/nf_mmseqs_trembl_query_dbs", f'{base_name}.mmseqs.trembl.in.pep')
#             with open(path_to_trembl_in_pep, 'r') as f:
#                 name_set = {_.rstrip()[1:] for _ in f if _.startswith('>')}
#             # now we can go through the dict and for any query in the name_set change the db value to trembl
#             # then write back out
#             new_df_dict = {}
#             count = 0
#             for k, v in dict_to_fix.items():
#                 if k in name_set:
#                     new_list = list(v)
#                     new_list[0] = 'trembl'
#                     new_df_dict[k] = new_list
#                     count += 1
#                 else:
#                     new_df_dict[k] = v
#             print(f'changed {count} out of {len(dict_to_fix.items())}')
#             pickle.dump(new_df_dict, open(path, 'wb'))

#     def _make_list_of_dict_paths(self):
#         cache_dir = "/home/humebc/projects/parky/strain_dnds/nf_go_annotations/cache"
#         return [os.path.join(cache_dir, _) for _ in os.listdir(cache_dir) if 'go_df_dict.p' in _]

if __name__ == "__main__":
    # CheekyFix()
    wotif = WriteOutPanzerInFasta()