import os
import sys
import pickle
from multiprocessing import Process, Queue
import subprocess

class RunCodeml:
    def __init__(self):
        self.species = sys.argv[1]
        self.summary_output_path = sys.argv[2]
        self.pickle_input_path = sys.argv[3]
        self.input_queue = Queue()

    def run_CODEML_analyses(self):
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

        num_proc = 1

        for i in range(num_proc):
            self.input_queue.put('STOP')

        list_of_processes = []
        for N in range(num_proc):
            p = Process(target=self.codeml_run_worker, args=(self.input_queue,))
            list_of_processes.append(p)
            p.start()

        for p in list_of_processes:
            p.join()

    def codeml_run_worker(self, input_queue):
        for dir_tup in iter(input_queue.get, 'STOP'):
            seq_file_path, ctrl_file_path, tree_file_path = dir_tup

            wkd = os.path.dirname(seq_file_path)
            # This may cause some problems
            # try to implement without the change of wkd.
            os.chdir(wkd)

            subprocess.run(['codeml', ctrl_file_path])