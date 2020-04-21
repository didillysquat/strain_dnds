#!/usr/bin/env python3

import sys

class WriteOutTrEMBLInFasta:
    def __init__(self):
        # The blast results as a dict
        # key will be fasta name line
        # Value will be the first result
        self.blast_out_dict = self._get_blast_out_list()
        # This will hold the contents of the new fasta
        self.new_pep_list = []
        # Read through the original .pep file and see which 
        # sequences didn't get a blast results
        # put these into the new pep
        self.all_seq_names = self._populate_new_pep_list()
        self.new_pep_file_name = sys.argv[1].replace('.sprot.out', '.trembl.in.pep')
        self._write_out_new_pep()
        

    def _get_blast_out_list(self):
        with open(sys.argv[1], 'r') as f:
            return [line.rstrip() for line in f]

    def _write_out_new_pep(self):
        with open(self.new_pep_file_name, 'w') as f:
            for line in self.new_pep_list:
                f.write(f'{line}\n')

    def self._populate_new_pep_list(self):
        """
        TODO work through the original pep file to see which seqs returned a blast result
        For those that didn't, add them to the new pep so that they can be blasted against
        the trembl db.
        """
        raise NotImplementedError
        