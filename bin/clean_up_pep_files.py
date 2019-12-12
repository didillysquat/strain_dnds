#!/usr/bin/env python3
"""
Here we will clean up the pep files so that they are exactly in fasta format
We will do this by removing everything after the name of the sequence
And by removing the asterix that are in the sequences
"""

import sys

class CleanPep:
    def __init__(self):
        self.path_to_files_to_clean = sys.argv[1:]

    def cleanup(self):
        # for each pep file, read it in, clean it and write it back out
        for pep_file in self.path_to_files_to_clean:
            with open(pep_file, 'r') as f:
                pep_list = [line.rstrip() for line in f]

            # create the clean pep file
            new_pep = []
            for line in pep_list:
                if line[0] == '>':
                    new_pep.append(line.split()[0])
                else:
                    new_pep.append(line.replace('*', ''))

            # Write out the new pep file
            with open(f'{pep_file}_clean', 'w') as f:
                for line in new_pep:
                    f.write(f'{line}\n')

cp = CleanPep()
cp.cleanup()