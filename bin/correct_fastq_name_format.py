#!/usr/bin/env python3
import sys
import os

class CorrectFastqNameFormat:
    def __init__(self):
        self.list_of_fastqs_to_correct = self._get_list_of_fastqs_to_correct()
        print(self.list_of_fastqs_to_correct)

    def _get_list_of_fastqs_to_correct(self):
        files = list(os.walk('.'))[0][2]
        return [file_name for file_name in files if '.fastq' in file_name]

    def correct_names(self):

        for fastq_path in self.list_of_fastqs_to_correct:
            print(f"Starting correction of {fastq_path}")
            new_fastq = []
            with open(fastq_path, 'r') as f:
                fastq_file = [line.rstrip() for line in f]

            for i in range(len(fastq_file)):
                current_line = fastq_file[i]
                if i%4 == 0 and current_line[0] == '@':
                    new_line = f"@{current_line.split(' ')[1]}"
                    new_fastq.append(new_line)
                    sys.stdout.write(f'\r{current_line} -> {new_line}')
                else:
                    new_fastq.append(current_line)

            if '_1.fastq' in fastq_path:
                out_path = fastq_path.replace('_1.fastq', '_nf_1.fastq')
            else:
                out_path = fastq_path.replace('_2.fastq', '_nf_2.fastq')
            with open(out_path, 'w') as f:
                for line in new_fastq:
                    f.write(f'{line}\n')

cfnf = CorrectFastqNameFormat()
cfnf.correct_names()