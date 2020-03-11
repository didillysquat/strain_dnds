#!/usr/bin/env python
"""
Writes out a cds alignment in phylip format, a control file, and the tree file.
Then runs the codeml analysis.
This script runs once per ortholog cds alignment.
"""
import os
import sys
import subprocess

class RunCODEML:
    def __init__(self):
        self.cds_alignment_path = sys.argv[1]
        self.orth_id = sys.argv[1].split('_')[0]
        self.tree_path_in = sys.argv[2]
        self.phylip_alignment = []
        self.fasta_file = []
        self.seq_file_path = f'{self.orth_id}_cds.phylip'
        self.out_file_path = f'{self.orth_id}_codeml_results.out'
        self.ctrl_path = f'{self.orth_id}_cds.ctrl'
        self.tree_path_out = f'{self.orth_id}_tree.nwk'


    def generate_block_phylip_alignments_for_CODEML(self):
        ''' This'''

        self.phylip_alignment = self._generate_phylip_from_fasta()
        self._write_out_phylip_file()
        self._make_and_write_control_file()
        self._write_out_tree_file()

    def _write_out_phylip_file(self):
        # write out the phylip file
        with open(self.seq_file_path, 'w') as f:
            for line in self.phylip_alignment:
                f.write(f'{line}\n')

    def _generate_phylip_from_fasta(self):
        temp_str = ''
        temp_list = []
        with open(self.cds_alignment_path, 'r') as f:
            fasta_file = [line.rstrip() for line in f]

        if len(fasta_file[1]) == 0:
            # if the fasta is empty then we need to log this outside
            raise RuntimeError('fasta empty')
        if len(fasta_file[1])%3 != 0:
            raise RuntimeError('not 3 modulus')
        else:
            for line in fasta_file:
                if line.startswith('>'):
                    temp_str = line[1:]
                else:
                    temp_list.append('{}    {}'.format(temp_str, line))
            # finally put in the header file
            temp_list.insert(0, f'\t{len(temp_list)} {len(fasta_file[1])}')

            return temp_list


    def _write_out_tree_file(self):
        with open(self.tree_path_in, 'r') as f:
            tree_file = [line.rstrip() for line in f]

        if len(tree_file) != 1:
            raise RuntimeError('Tree file is not len 1')

        with open(self.tree_path_out, 'w') as f:
            f.write(tree_file[0])

    def _make_and_write_control_file(self):
        ctrl_file = [
        f'seqfile = {self.seq_file_path}',
        f'treefile = {self.tree_path_out}',
        f'outfile = {self.out_file_path}',
        'runmode = -2',
        'seqtype = 1',
        'CodonFreq = 2',
        'ndata = 1',
        'clock = 0',
        'model = 0',
        'NSsites = 0',
        'icode = 0',
        'fix_omega = 0',
        'omega = .4',
        'cleandata = 0',
        'verbose = 1'
        ]

        with open(self.ctrl_path, 'w') as f:
            for line in ctrl_file:
                f.write('{}\n'.format(line))

    def run_codeml(self):
        subprocess.run(args=['codeml', self.ctrl_path])

cml = RunCODEML()
# We will write out a status file whether this was successful or not
# This is so that the process will be forced to run as the .out file is optional
# We will make it so that the <id>_status.txt id mandatory 
try:
    cml.generate_block_phylip_alignments_for_CODEML()
    cml.run_codeml()
    with open(f'{cml.orth_id}_status.txt', 'w') as f:
        f.write('0\n')
except:
    # We have set the script to raise RuntimeError if the cds alignment is
    # empty or is not divisible by 3
    # If this happens we want to exit out
    # There will not have been a codeml output file produced
    with open(f'{cml.orth_id}_status.txt', 'w') as f:
        f.write('1\n')
    pass
