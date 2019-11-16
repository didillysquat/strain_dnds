# output argument will always be sonicparanoid/output
# input will always be sonicparanoid

import os
import sys
import shutil
import subprocess

class ParanoidAnalysis:
    def __init__(self):
        self.threads = int(sys.argv[1])
        self.input_dir = os.path.abspath("sonicparanoid")
        self.output_dir = os.path.join(self.input_dir, "output")
        self.runs_dir = os.path.join(self.output_dir, "runs")
        self.parkinson_dir = os.path.join(self.output_dir, "runs", "parkinson")

    def run(self):
        self._del_previous_runs()
        self._run_sonicparanoid()
        self._rename_run()

    def _del_previous_runs(self):
        print("checking to see if sonicparanoid/output/runs/parkinson/ exists")
        if os.path.exists(self.parkinson_dir):
            print(f"{self.parkinson_dir} already exists")
            print('Deleting')
            shutil.rmtree(self.parkinson_dir)
            print("Done")
        else:
            print("sonicparanoid/output/runs/parkinson/ does not exist")

    def _run_sonicparanoid(self):
        print(f"Running sonicparanoid analysis with {self.threads} threads")
        print(f"sonicparanoid -i {self.input_dir} -o {self.output_dir} -t {self.threads}")
        subprocess.run(['sonicparanoid', '-i', self.input_dir, '-o', self.output_dir, '-t', self.threads])

    def _rename_run(self):
        """
        because we didn't specify the run name (this was causing problems) we will rename the output
        that is automatically by sonicparanoid from the current date and time. We will call it parkinson
        so that we know what the output is for the snakemake file. We could use a config file at a later date
        to specify what this will be called both here and in the snakemakefile
        """
        # There should be only one directlry in the runs directory
        # check to see if there are two and raise an excpetion is there are
        dir_contents = list(os.walk(self.runs_dir))[0][1]
        num_dirs = len(dir_contents)
        if num_dirs == 0:
            raise RuntimeError(f"There is no directory in {self.runs_dir}")
        elif num_dirs > 1:
            raise RuntimeError(f"There are {num_dirs} directories in {self.runs_dir}")
        else:
            # There is only one dir that needs to be renamed
            src = os.path.join(self.runs_dir, dir_contents[0])
            os.rename(src=src, dst=self.parkinson_dir)

pa = ParanoidAnalysis()
pa.run()