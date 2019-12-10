#!/usr/bin/env python
import sys
import subprocess
import ntpath

class RunModelTest:
    def __init__(self):
        self.cropped_aligned_fasta = sys.argv[1]
        self.orth_id = ntpath.basename(sys.argv[1]).split('_')[0]
        
    def run_model_test(self):
        # first check to see if there is actually anything to run the test on
        with open(self.cropped_aligned_fasta, 'r') as f:
            fasta_file = [line.rstrip() for line in f]
        # write out a empty .out file and return
        if not fasta_file[1]:
            self._write_out_blank_test_out()
            return
        
        # otherwise, run the prottest
        subprocess.run(['modeltest-ng', '-i', self.cropped_aligned_fasta, '-d', 'aa', '-p', '1', '-o', self.cropped_aligned_fasta.replace('_cropped_aligned_aa.fasta', '_prottest_result'), '-m', '+LG4M,+LG4X,+GTR'])

    def _write_out_blank_test_out(self):
        with open(self.cropped_aligned_fasta.replace('_cropped_aligned_aa.fasta', '_prottest_result.out'), 'w') as f:
                f.write('')

rmt = RunModelTest()
rmt.run_model_test()