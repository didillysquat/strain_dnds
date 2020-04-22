#!/usr/bin/env python3

import sys
import os
from collections import defaultdict

class WriteOutTrEMBLInFasta:
    def __init__(self):
        # A set of the seq names that returned a result
        self.blast_out_name_set = self._make_blast_out_name_set()
        # List of the full paths in which to search for the pep file of interest
        self.pep_input_directories = sys.argv[2].split(',')
        self.name_of_pep_file = sys.argv[1].replace('.sprot.out', '.pep')
        self.pep_file_dict = None
        for dir in self.pep_input_directories:
            if os.path.isfile(os.path.join(dir, self.name_of_pep_file)):
                found = True
                with open(os.path.join(dir, self.name_of_pep_file), 'r') as f:
                    pep_file_list = [line.rstrip() for line in f]
                    self.pep_file_dict = {pep_file_list[i].split()[0][1:]:pep_file_list[i+1] for i in range(0, len(pep_file_list), 2)}
                    # self.pep_file_dict = {f[i]:f[i+1], for i in range(0, len(f), 2)}
                break
        if self.pep_file_dict is None:
            raise RuntimeError("Couldn't locate original .pep file")
        
        # Read through the original .pep file and see which 
        # sequences didn't get a blast results
        # put these into the new pep
        self.new_pep_list = self._make_new_pep_list()
        self.new_pep_file_name = sys.argv[1].replace('.sprot.out', '.trembl.in.pep')
        self._write_out_new_pep()
        

    def _make_blast_out_name_set(self):
        with open(sys.argv[1], 'r') as f:
            blast_out_list = [line.rstrip() for line in f]
        
        blast_out_name_set = set()
        for line in blast_out_list:
            blast_out_name_set.add(line.split('\t')[0])
        return blast_out_name_set

    def _write_out_new_pep(self):
        with open(self.new_pep_file_name, 'w') as f:
            for line in self.new_pep_list:
                f.write(f'{line}\n')
        

    def _make_new_pep_list(self):
        """
        TODO work through the original pep file to see which seqs returned a blast result
        For those that didn't, add them to the new pep so that they can be blasted against
        the trembl db.
        """
        new_pep_list = []
        for k, v in self.pep_file_dict.items():
            if k not in self.blast_out_name_set:
                # Then we had no blast result for this seq and it needs to be added to the
                # new pep file
                new_pep_list.extend([f'>{k}', v])
            else:
                pass
        return new_pep_list

if __name__ == "__main__":
    wotif = WriteOutTrEMBLInFasta()