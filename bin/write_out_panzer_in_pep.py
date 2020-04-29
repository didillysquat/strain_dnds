#!/usr/bin/env python3
# This script will take in the output mmseq files from running the nr search
# and look to see which of the matches allow us to associate a GO.
# If we are not able to associate a GO then we will write out the query sequence
# into a new pep that will be run against the panzer server.
# For each query that we do get a hit for, we will add it to the go_df_dict
# that was created and pickled out as part of the processing of the sprot out
# results.
# To do the mapping we will first attempt to map to an ncbi gene id
# using gene2annotation. We will then try to map this geneid
# to a GO using gene2go.
# We are expecting the return from this to be very low.
# We should also convert the go_df_dict into a final dataframe that can be output
# for john.

import sys
import os
from collections import defaultdict
import pandas as pd
import pickle

class WriteOutTrEMBLInFasta:
    def __init__(self):
        # A set of the seq names that returned a result
        # And the actual output match file as a list
        self.mmseqs_out_name_set, self.mmseqs_match_list = self._make_mmseqs_out_name_set()
        # List of the full paths in which to search for the pep file of interest
        self.pep_input_directory = sys.argv[2]
        self.original_pep_file_name = sys.argv[3].split('/')[-1]
        self.cache_dir_base = sys.argv[4]
        self.base_name = self.original_pep_file_name.replace('.mmseqs.nr.in.pep', '')
        self.uniprot_to_go_dict_cache_path = os.path.join(self.cache_dir_base, f'{self.base_name}_uniprot_to_go_dict.p')
        self.uniprot_to_go_dict = pickle.load(open(self.uniprot_to_go_dict_cache_path, 'rb'))
        self.pep_file_dict = None
        self.df_columns = []
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
        self.gene2accession_path = "/share/databases/gene2accession/gene2accession"
        

        # now we need to go through the mmseqs file and see which matches can be associated to a GO
        self._log_go_associations()

        self.new_pep_file_name = self.original_pep_file_name.replace('.mmseqs.nr.in.pep', '.mmseqs.panzer.in.pep')
        self._write_out_new_pep()
        # Pickle out the go_df_dict
        pickle.dump(self.go_df_dict, open(self.go_df_dict_path, 'wb'))
        
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
            # If there are less than 2000 match lines left to process, then let's just
            # request all the matches at once, to save on time.
            if len(sub_match_query_set) + len(sub_match_list_for_next_iteration) < 20000:
                match_accessions = []
                match_accessions.extend([_.split('\t')[1] for _ in sub_match_list])
                match_accessions.extend([_.split('\t')[1] for _ in sub_match_list_for_next_iteration])
                match_accessions = set(match_accessions)
            else:
                match_accessions = set([_.split('\t')[1] for _ in sub_match_list])
            accessions_to_retrieve = [_ for _ in match_accessions if _ not in self.uniprot_to_go_dict.keys()]
            # Then we will work in chunks to reduce the chances of a connection reset etc.
            if accessions_to_retrieve:
                self._update_uniprot_to_go_dict(accessions_to_retrieve)
            
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
                    go_list = self.uniprot_to_go_dict[match]
                    if go_list != 'no_annotation':
                        self.go_df_dict[line.split('\t')[0]] = ['sprot', match, line.split('\t')[2], match_number, go_list]
                except KeyError:
                    raise RuntimeError(f'{match} was not found in the uniprot_to_go_dict')

    def _update_uniprot_to_go_dict(self, accessions_to_retrieve):
        self.input_queue = Queue()
        self.output_queue = Queue()
        self.num_proc = 10
        for chunk in self._chunks(accessions_to_retrieve, 250):
            self.input_queue.put(chunk)
        
        for n in range(self.num_proc):
            self.input_queue.put('STOP')

        all_processes = []
        for n in range(self.num_proc):
            p = Process(target=self._get_request_results, args=())
            all_processes.append(p)
            p.start()
        
        done_count = 0
        while done_count < self.num_proc:
            result_list = self.output_queue.get()
            if result_list == 'DONE':
                done_count += 1
            else:
                for line in result_list[1:]:
                    go_list = line.split('\t')[0]
                    if not go_list:
                        self.uniprot_to_go_dict[line.split('\t')[1]] = 'no_annotation'
                    else:
                        self.uniprot_to_go_dict[line.split('\t')[1]] = ''.join(go_list.split(' '))
            # Now pickle out the dict
            print(f'Main thread: successfully added {len(result_list) - 1} items to the uniprot_to_go_dict, total now {len(self.uniprot_to_go_dict.items())}')
            pickle.dump(self.uniprot_to_go_dict, open(self.uniprot_to_go_dict_cache_path, 'wb'))

        for p in all_processes:
            p.join()
        foo = 'ar'

    def _get_request_results(self):
        for chunk in iter(self.input_queue.get, 'STOP'):
            accessions_to_retrieve = ' '.join(chunk)
            response = self._make_request(accessions_to_retrieve)
            go_info = StringIO(response)
            self.output_queue.put([_.rstrip() for _ in go_info.readlines()])
        self.output_queue.put('DONE')


    def _chunks(self, lst, n):
        """Yield successive n-sized chunks from lst."""
        for i in range(0, len(lst), n):
            yield lst[i:i + n]
    
    def _make_request(self, accessions_to_retrieve):
        while True:
            url = 'https://www.uniprot.org/uploadlists/'
            params = {
            'from': 'ACC+ID',
            'to': 'ACC',
            'format': 'tab',
            'query': accessions_to_retrieve,
            'columns': 'go-id'
            }
            data = urllib.parse.urlencode(params)
            data = data.encode('utf-8')
            headers = {'User-Agent': 'Benjamin Hume', 'From': 'benjamin.hume@kaust.edu.sa'}
            # req = urllib.request.Request(url, data, headers)
            # try:
            print(f'{current_process().name}: attempting request')
            try:
                r = requests.get(url=url, params=params, headers=headers)
                if r.status_code != 200:
                    print(f"{current_process().name}: status code {r.status_code} returned")
                    print('Waiting 30s before retrying connection')
                    time.sleep(30)
                    continue
                print(f"{current_process().name}: request successful!")
                return r.text
            except requests.exceptions.ConnectionError as e:
                print(f'{current_process().name}: {e}')
                print(f'{current_process().name}: Waiting 30s before retrying connection')
                time.sleep(30)
                continue
            

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

if __name__ == "__main__":
    wotif = WriteOutTrEMBLInFasta()