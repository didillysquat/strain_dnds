import os
import sys
import pickle
from multiprocessing import Process, Queue
import subprocess

class RunCodeml:
    def __init__(self):
        self.species = sys.argv[1]
        self.pickle_input_path = sys.argv[2]
        self.num_proc = 12
        self.input_queue = Queue()
        self.summary_out_path = f"/home/humebc/projects/parky/breviolum_transcriptomes/codeml/" \
                                f"{self.species}/codeml_run_summary.txt"

    def run_codeml_analyses(self):
        ''' We should read in the tuples that contain the seq files and cntrl file
        from the pickled files that were written out from generate_master_phlip_alignments_for_CODEML
        We should then start an MP list with each of these tuples in
        In the worker we should go to each directory and start an anlysis
        The 20000 sequences were broken into 20 chunks so we should start 20 subprocess instances and
        have a one CODEML analysis run in each'''

        tup_of_dirs_list = pickle.load(
            open(self.pickle_input_path, "rb"))

        for tup in tup_of_dirs_list:
            self.input_queue.put(tup)

        for i in range(self.num_proc):
            self.input_queue.put('STOP')

        list_of_processes = []
        for N in range(self.num_proc):
            p = Process(target=self.codeml_run_worker, args=(self.input_queue,))
            list_of_processes.append(p)
            p.start()

        for p in list_of_processes:
            p.join()

        with open(self.summary_out_path, 'w') as f:
            f.write('DONE')

    def codeml_run_worker(self, input_queue):
        for dir_tup in iter(input_queue.get, 'STOP'):
            seq_file_path, ctrl_file_path, tree_file_path = dir_tup

            wkd = os.path.dirname(seq_file_path)
            # This may cause some problems
            # try to implement without the change of wkd.
            os.chdir(wkd)
            print('Starting codeml')
            print(ctrl_file_path)
            subprocess.run(args=['/usr/bin/codeml', ctrl_file_path])

rc = RunCodeml()
rc.run_codeml_analyses()