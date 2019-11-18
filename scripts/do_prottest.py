import os
from multiprocessing import Queue, Process
import subprocess
import sys


class DoProttest:
    def __init__(self):
        self.species = sys.argv[1]
        self.threads = int(sys.argv[2])
        self.base_input_dir = os.path.join('/home/humebc/projects/parky/breviolum_transcriptomes/local_alignments', self.species)
        self.prottest_path = '/home/humebc/phylogeneticSoftware/prottest/prottest-3.4.2/prottest-3.4.2.jar'
        self.list_of_aa_alignments_mp_queue = Queue()
        self._populate_list_of_aa_alignments()
        self.alignment_not_available_list = []

    def _populate_list_of_aa_alignments(self):
        list_of_orth_group_dirs = list(os.walk(self.base_input_dir))[0][1]
        print("Populating the list of aa alignments mp queue")
        for group_dir_name in list_of_orth_group_dirs:
            sys.stdout.write(f"\r{group_dir_name}")
            path_to_aa_alginment = os.path.join(self.base_input_dir, group_dir_name, f'{group_dir_name}.cropped_aligned_aa.fasta')
            if os.path.exists(path_to_aa_alginment):
                self.list_of_aa_alignments_mp_queue.put(path_to_aa_alginment)
            else:
                self.alignment_not_available_list.append(path_to_aa_alginment)

    def do_mp_prottest(self):
        # put in the STOPs
        for _ in range(self.threads):
            self.list_of_aa_alignments_mp_queue.put('STOP')

        all_processes = []

        # Then start the workers
        for _ in range(self.threads):
            p = Process(target=self._prottest_worker, args=(self.list_of_aa_alignments_mp_queue,))
            all_processes.append(p)
            p.start()

        for p in all_processes:
            p.join()

        self._write_prottest_summary()

    def _write_prottest_summary(self):
        with open(os.path.join(self.base_input_dir, 'protein_models_summary.txt'), 'w') as f:
            f.write(f'{len(self.alignment_not_available_list)} alignments were not found in their respective directories')

    def _prottest_worker(self, input_queue):
        for aa_alignment_path in iter(input_queue.get, 'STOP'):
            output_path = aa_alignment_path.replace('.cropped_aligned_aa.fasta', '_prottest_result.out')
            # if the test has already been done
            if os.path.isfile(output_path):
                print(f'Already exists {output_path}')
                continue
            sys.stdout.write(f'\rRunning prottest: {aa_alignment_path}')
            # perform the prottest
            subprocess.run(
                ['java', '-jar', self.prottest_path, '-i', aa_alignment_path, '-o', output_path, '-all-distributions', '-all'])

dp = DoProttest()
dp.do_mp_prottest()