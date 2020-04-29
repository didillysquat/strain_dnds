"""
Purpose of this script is to create a pickled out dictionary
that maps all of matches from the to a gene ID if there is one
available using the gene2accession.
We have split the gene2accession file into 200 files.
We will create 200 processes and in each one we will
read in one of the files and create a dict of nr accesion
to gene id. In to each of these process we will pass the master
set that contains all of the nr accessions that were matches
to the gene ids from the mmseqs nr out.
Then in each process we will perform an instersect of the two sets,
subsample both the dicts accordingly and finally make a mapping
between the two dicts.
We will collect all of the dicts outside of the MP framework
before combining them and pickle ing out as a single dict
that we can then use in the follow up processes.
"""

import sys
import os
from multiprocessing import Queue, Process

class MakeGene2AccessionMapping:
    def __init__(self):
        self.chunk_dir = sys.argv[1]
        self.chunk_paths = [os.path.join(self.chunk_dir, _) for _ in os.listdir(self.chunk_dir)]
        self.nr_out_results_dir = sys.argv[2]
        self.nr_out_paths = [os.path.join(self.nr_out_results_dir, _) for _ in os.listdir(self.nr_out_results)]
        self.master_match_set = self._make_master_match_set()
        self.input_queue = Queue()
        self.num_proc = 200
        self._setup_and_run_mp()

    def _setup_and_run_mp(self):
        for chunk_path in self.chunk_paths:
            self.input_queue.put(chunk_path)
        for n in range(self.num_proc):
            self.input_queue.put('STOP')
        
        all_processes = []
        for n in range(self.num_proc):
            p = Process(target=self._get_request_results, args=(self.master_match_set))
            all_processes.append(p)
            p.start()

        done_count = 0
        master_map = {}
        while done_count < self.num_proc:
            sub_map = self.output_queue.get()
            if sub_map == 'DONE':
                done_count += 1
            else:
                master_map.update(sub_map)

        for p in all_processes:
            p.join()

        # Now we can pickle out the dict
        pickle.dump(master_map, open('gene2accession.p', 'wb'))

    def _make_sub_mapping(self, match_set):
        for chunk_path in iter(self.input_queue.get, 'STOP'):
            # Make the raw mapping dict
            with open(chunk_path, 'r') as f:
                map_dict = {_.split('\t')[1]: _.split('\t')[5] for _ in f}
            # Get the intersection
            intersect = set(map_dict.keys()).intersection(match_set)
            # slim down match_set and map_dict
            return {k:v for k, v in map_dict.items() if k in intersect}
        self.output_queue.put('DONE')
            
    def _make_master_match_set(self):
        master_set = set()
        for out_path in self.nr_out_paths:
            with open(output_path, 'r') as f:
                lines = [_.rstrip() for _ in f]
                master_set.update([_.split('\t')[1] for _ in lines])
        return master_set

if __name__ == "__main__":
    MakeGene2AccessionMapping()