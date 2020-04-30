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
from multiprocessing import Queue, Process, current_process
import pickle
from collections import defaultdict

class MakeGene2AccessionMapping:
    def __init__(self):
        self.chunk_dir = sys.argv[1]
        self.chunk_paths = [os.path.join(self.chunk_dir, _) for _ in os.listdir(self.chunk_dir)]
        self.nr_out_paths = [_ for _ in os.listdir('.') if 'blast_format' in _]
        self.master_match_set = self._make_master_match_set()
        self.input_queue = Queue()
        self.output_queue = Queue()
        self.num_proc = 50
        self._setup_and_run_mp()

    def _setup_and_run_mp(self):
        for chunk_path in self.chunk_paths:
            self.input_queue.put(chunk_path)
        for n in range(self.num_proc):
            self.input_queue.put('STOP')
        
        all_processes = []
        for n in range(self.num_proc):
            p = Process(target=self._make_sub_mapping, args=(self.master_match_set,))
            all_processes.append(p)
            p.start()


        # want to take into account that a given 
        done_count = 0
        master_map = defaultdict(set)
        while done_count < self.num_proc:
            sub_map = self.output_queue.get()
            if sub_map == 'DONE':
                done_count += 1
            else:
                # for k, v in sub_map.items():
                #     master_map[k].update
                for k in sub_map.keys():
                    if k in master_map:
                        print(f'key is {k}')
                        print(f'master_map[k] was: {master_map[k]}')
                        master_map[k].update(sub_map[k])
                        print(f'master_map[k] now: {master_map[k]}')
                    else:
                        master_map[k].update(sub_map[k])
                # master_map.update(sub_map)

        for p in all_processes:
            p.join()

        # Now we can pickle out the dict
        pickle.dump(master_map, open('gene2accession.p', 'wb'))

    def _make_sub_mapping(self, match_set):
        chunk_num = 0
        for chunk_path in iter(self.input_queue.get, 'STOP'):
            chunk_num += 1
            print(f'{current_process().name}: chunk num {chunk_num}')
            # Make the raw mapping dict
            # We will take into account that a given accession may be linked to several gene ids
            map_dict = defaultdict(set)
            with open(chunk_path, 'r') as f:
                for line in f:
                    if line.split('\t')[5] != '-':
                        map_dict[line.split('\t')[5]].add(line.split('\t')[1])
            # Get the intersection
            intersect = set(map_dict.keys()).intersection(match_set)
            # slim down match_set and map_dict
            self.output_queue.put({k:v for k, v in map_dict.items() if k in intersect})
        self.output_queue.put('DONE')
            
    def _make_master_match_set(self):
        master_set = set()
        for out_path in self.nr_out_paths:
            with open(out_path, 'r') as f:
                lines = [_.rstrip() for _ in f]
                master_set.update([_.split('\t')[1] for _ in lines])
        return master_set

if __name__ == "__main__":
    MakeGene2AccessionMapping()