#!/usr/bin/env python3
# This script will take in the output mmseq files from running the sprot search
# and look to see which of the matches allow us to associate a GO.
# If we are not able to associate a GO then we will write out the query sequence
# into a new pep that will be run against the trembl database.
# For each query that we do get a hit for, we will record add a row to the df
# The df will have columns:
# db, query_name, match_accession, match_e_value, best_hit_number, GO_accessions
# db=the database that we got a hit from
# qyery_name = the name of the query sequence
# match_accession = the matching accession that we got a GO annotation from
# match_e_value = the e value of the accession that matched the sequence
# best_hit number will be how far down the list of matches we found the GO_accession (1=first match)
# GO_accession will be a comma seperated list of the GO terms associated

import sys
import os
from collections import defaultdict
import pandas as pd
import urllib.parse
import urllib.request
import pickle
from io import StringIO

class WriteOutTrEMBLInFasta:
    def __init__(self):
        # A set of the seq names that returned a result
        # And the actual output match file as a list
        self.mmseqs_out_name_set, self.mmseqs_match_list = self._make_mmseqs_out_name_set()
        # List of the full paths in which to search for the pep file of interest
        self.pep_input_directories = sys.argv[2].split(',')
        self.name_of_pep_file = sys.argv[3].split('/')[-1]
        self.cache_dir_base = sys.argv[4]
        self.base_name = self.name_of_pep_file.replace('_longest_iso_orfs.single_orf.pep', '')
        self.uniprot_to_go_dict_cache_path = os.path.join(self.cache_dir_base, f'{self.base_name}_uniprot_to_go_dict.p')
        self.uniprot_to_go_dict = defaultdict(list)
        self.pep_file_dict = None
        self.go_associations_dict = defaultdict(list)
        for dir in self.pep_input_directories:
            if os.path.isfile(os.path.join(dir, self.name_of_pep_file)):
                found = True
                with open(os.path.join(dir, self.name_of_pep_file), 'r') as f:
                    pep_file_list = [line.rstrip() for line in f]
                    self.pep_file_dict = {pep_file_list[i].split()[0][1:]:pep_file_list[i+1] for i in range(0, len(pep_file_list), 2)}
                break
        if self.pep_file_dict is None:
            raise RuntimeError("Couldn't locate original .pep file")
        
        # the df that will hold the GO annotations that we were able to find.
        # We will make this as a dict to start with and then create from there at the end
        self._go_df_dict = {}
        self._go_df_headers = ['db', 'match', 'e_value', 'hit_number', 'GO_accessions']
        # Read through the original .pep file and see which 
        # sequences didn't get a blast results
        # put these into the new pep
        self.new_pep_list = self._add_no_match_seqs_to_pep()
        
        # now we need to go through the mmseqs file and see which matches can be associated to a GO
        self._log_go_associations()

        self.new_pep_file_name = self.name_of_pep_file.replace('_longest_iso_orfs.single_orf.pep', '.mmseqs.trembl.in.pep')
        self._write_out_new_pep()
        
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
                if query in self._go_df_dict:
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
            match_accessions = set([_.split('\t')[1] for _ in sub_match_list])
            self._update_uniprot_to_go_dict(match_accessions)
            # Now we have the mappings that we need to see if there are go associations
            for line in self.sub_match_list:
                match = line.split('\t')[1]
                go_list = self.uniprot_to_go_dict[match]
                if go_list:
                    # log the association
                    self.go_associations_dict[line.split('\t')[0]] = ['sprot', match, line.split('\t')[2], match_number, go_list]
                # nothing to do if no succesful association
            if not sub_match_list_for_next_iteration:
                break
            match_number += 1


    def _update_uniprot_to_go_dict(self, set_of_uniprot_accesions):
        """
        This function will take in a request for a set_of_uniprot_accessions
        It will load up our locally cached dictionary of uniprot accessions to GO
        terms and it will see which of the accesions in the input set
        are not already contained in our local dictionary.
        If there are some that are not contained, then it will do a request for them,
        update the dictionary and pickle out.
        We will keep a cache per transciptome as we may hit issues with concurrent access
        otherwise.
        """
        # First check to see if the cache already exists
        if os.path.exists(self.uniprot_to_go_dict_cache_path):
            self.uniprot_to_go_dict = pickle.load(open(self.uniprot_to_go_dict_cache_path, 'rb'))
        # else we have already initiated it as a default dict
        #get a set of the accessions that we need to get
        accessions_to_retrieve = ' '.join([_ for _ in set_of_uniprot_accesions if _ not in self.uniprot_to_go_dict.keys()][:10])
        if accessions_to_retrieve:
            #TODO if there are some to get then get them, else just return the dict as it is.
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
            req = urllib.request.Request(url, data, headers)
            with urllib.request.urlopen(req) as f:
                response = f.read()
            for line for 
            go_info = StringIO(response.decode('utf-8')))
            for line in go_info[1:]:
                go_list = line.split('\t')[0]
                if not go_list:
                    continue
                else:
                    self.uniprot_to_go_dict[line.split('\t')[1]] = ''.join(go_info.split(' '))
                    # TODO we are here.
            foo = response.decode('utf-8')
            print(foo)
            foo = 'bar'
            #TODO process the output so that we add to the dictionary, pickle out and
            # update self.uniprot_to_go_dict
            # return 
        else:
            # there are no new uniprot seqs to get the go terms for and we
            # can work with the dict that we already have
            return 

    def _make_goa_uniprot_all_df(self):
        # We 
        # The .gaf file is massive so lets try to be smart about finding matches
        # first collect all a big set of all of the match accessions
        first_10_matches = ' '.join(list(set([_.split('\t')[1] for _ in self.mmseqs_match_list[:10]])))
        

        # url = "http://www.uniprot.org/uniprot/?query=Q9SLA1,Q3EAF9,P40371,Q9SD02,Q69VD9,P35182,Q8N819,Q6AUQ4,Q9LNF4,Q9FYN7&format=tab&columns=id%2Cgo"
        url = 'https://www.uniprot.org/uploadlists/'

        params = {
        'from': 'ACC+ID',
        'to': 'ACC',
        'format': 'tab',
        'query': first_10_matches,
        'columns': 'go-id'
        }

        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        headers = {'User-Agent': 'Benjamin Hume', 'From': 'benjamin.hume@kaust.edu.sa'}
        req = urllib.request.Request(url, data, headers)
        with urllib.request.urlopen(req) as f:
            response = f.read()
        print(response.decode('utf-8'))

        
        match_set = set([_.split('\t')[1] for _ in self.mmseqs_match_list])
        
        # with open(self.path_to_goa_uniprot_all_gaf, 'r') as f:
        #     for line in f:
        #         if line.startswith('!'):
        #             continue
        #         else:


        df = pd.read_table(self.path_to_goa_uniprot_all_gaf, skiprows=[0,1,2,3,4,5,6,7])
        foo = 'bar'

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