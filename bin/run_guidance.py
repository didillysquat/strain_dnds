#!/usr/bin/env python
import subprocess
import sys

class RunGuidance:
    def __init__(self):
        self.guidance_perl_script_full_path = self._init_guidance_perl_script_full_path()
        self.input_fasta = sys.argv[1]
        self.orth_group_id = self.input_fasta.split('_')[0]
        self.guidance_run_attempts = 0

    def _init_guidance_perl_script_full_path(self):
        proc = subprocess.run(['which', 'guidance'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        guidance_path = proc.stdout.decode("utf-8").rstrip()
        with open(guidance_path, 'r') as f:
            for line in f:
                if 'guidance.pl' in line:
                    return line.split(' ')[1]
    
    def run_guidance(self):
        try:
            subprocess.run(
                        ['perl', self.guidance_perl_script_full_path, '--seqFile', self.input_fasta, '--msaProgram', 'MAFFT',
                        '--seqType', 'codon', '--outDir', '.', '--bootstraps', '10',
                        '--outOrder', 'as_input', '--colCutoff', '0.6', '--dataset', self.orth_group_id, '--proc_num', '1'])
        except:
            if self.guidance_run_attempts < 5:
                self.guidance_run_attempts += 1
                self.run_guidance()
            else:
                raise RuntimeError(f'Despite trying to rerun the Guidance analysis 5 times, there seems to be something wrong with {self.input_fasta}')

rg = RunGuidance()
rg.run_guidance()