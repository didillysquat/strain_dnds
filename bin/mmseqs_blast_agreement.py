"""
The purpose of this script is to see if the results from the mmseqs and blast give good agreement
For each sample we will read in the blast output and create a dict of query to result for the blast.
We will then compare this to the output of the mmseqs.
"""

import os
from collections import Counter

class MMSEQ_BLAST_Compare:
    def __init__(self):
        self.blast_out_dir = "/home/humebc/projects/parky/strain_dnds/nf_blast_results"
        self.mmseqs_out_dir = "/home/humebc/projects/parky/strain_dnds/nf_mmseqs_sprot_results"
        self.list_of_blast_out_files = self._get_list_of_blast_out_files()
        self.list_of_mmseq_out_files = self._get_list_of_mmseq_out_files()
        self.m_dict = {}
        self.current_query = None
        self.found_match = False
        self.match_count = 1
        self.blast_result = None
        self.skipping = False
        self.mmseq_set = set()
        self.r_dict = {}

    def compare(self):
        top_5_list = []
        top_1_list = []
        for blast_out_path in self.list_of_blast_out_files:
            self.m_dict = {}
            self.current_query = None
            self.found_match = False
            self.match_count = 1
            self.blast_result = None
            self.skipping = False
            self.mmseq_set = set()
            self.r_dict = {}
            with open(blast_out_path, 'r') as f:
                blast_list = [_.rstrip() for _ in f]

            # First log the blast output into a dict
            
            last_query = None
            for line in blast_list:
                query = line.split()[0]
                if query == last_query:
                    # Then we've already processed this and we just pass on
                    continue
                else:
                    # Then this is a best match of the query and we log
                    self.r_dict[query] = line.split()[2].split('|')[1]
                    last_query = query

            # Then check this dict against the mmseq results
            
            if 'b_min' in blast_out_path:
                mmseq_file_name = "SRR_b_min_c.mmseqs.out.blast_format"
            elif 'b_psyg' in blast_out_path:
                mmseq_file_name = "SRR_b_psyg_c.mmseqs.out.blast_format"
            elif 'BreviolumB5' in blast_out_path:
                mmseq_file_name = "BreviolumB5_Sradians.mmseqs.out.blast_format"
            elif 'Breviolumfaviinorum' in blast_out_path:
                mmseq_file_name = "Breviolumfaviinorum_Pclivosa.mmseqs.out.blast_format"
            else:
                mmseq_file_name = blast_out_path.split('/')[-1].split('_')[0] + ".mmseqs.out.blast_format"
            mmseq_out_path = os.path.join(self.mmseqs_out_dir, mmseq_file_name)
            with open(mmseq_out_path, 'r') as f:
                mmseq_list = [_.rstrip() for _ in f]

            # Go through line and check which result matches the one that matched the dict
            self.mmseq_set = set()
            for line in mmseq_list:
                query = line.split()[0]
                self.mmseq_set.add(query)
                if query == self.current_query and not self.skipping:
                    self._log_match(line)
                elif query != self.current_query:
                    # Then we need to reset
                    self.skipping = False
                    self.match_count = 1
                    self.current_query = query
                    self.found_match = False
                    if query in self.r_dict:
                        self.blast_result = self.r_dict[query]
                        self._log_match(line)
                    else:
                        # Then there was an mmseq result but no blast result.
                        # We need to just continue to the next query.
                        self.current_query = query
                        self.skipping = True
                        continue  
                else:
                    continue
            
            # At this point we can report for a given transcript
            print(f"{mmseq_file_name.split('.')[0]}:")
            print(f"\t{len(self.r_dict.keys())} blast results")
            print(f"\t{len(self.mmseq_set)} mmseq results")
            print(f"\t{len(self.m_dict.keys())} of the blast results had mmseq results")
            print(f"\t{len([_  for _ in self.mmseq_set if _ not in self.m_dict])} of mmseq results were not found in the blast")
            values_set = sorted(list(set(self.m_dict.values())))
            counter = Counter(self.m_dict.values())
            top_5 = 0
            for val in values_set[:5]:
                print(f"\t{val}: found {counter[val]} times ({(counter[val]/len(self.m_dict.values())) * 100}%)")
                top_5 += (counter[val]/len(self.m_dict.values())) * 100
            print(f"\ttop 5 represent: {top_5}")
            print('\n')
            top_5_list.append(top_5)
            top_1_list.append((counter[1]/len(self.m_dict.values())) * 100)
        # At this point we can plot up or something the m_dict
        print(f"Finally, for the 6.5 s value on average the top5 matched {sum(top_5_list)/len(top_5_list)}.")
        print(f"Finally, for the 6.5 s value on average the top1 matched {sum(top_1_list)/len(top_1_list)}.")

    
    def _log_match(self, line):
        if not self.found_match:
            if line.split()[1] == self.blast_result:
                self.m_dict[self.current_query] = self.match_count
                self.found_match = True
            else:
                self.match_count += 1

    def _get_list_of_blast_out_files(self):
        blast_out_path_list = []
        for file_path in os.listdir(self.blast_out_dir):
            blast_out_path_list.append(os.path.join(self.blast_out_dir, file_path))
        return blast_out_path_list

    def _get_list_of_mmseq_out_files(self):
        mmseq_out_path_list = []
        for file_path in os.listdir(self.mmseqs_out_dir):
            mmseq_out_path_list.append(os.path.join(self.mmseqs_out_dir, file_path))
        return mmseq_out_path_list

MMSEQ_BLAST_Compare().compare()