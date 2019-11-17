import os
import sys
from multiprocessing import Queue, Process
import ntpath
import subprocess
import pandas as pd
import shutil

class GAlign:
    def __init__(self):
        self.species = sys.argv[1]
        print(f'argvs 1 = {sys.argv[1]}')
        self.base_alignment_dir = os.path.join('/home/humebc/projects/parky/breviolum_transcriptomes/local_alignments', self.species)
        self.threads = 10
        # self.threads = int(sys.argv[3])
        self.mp_queue = Queue()
        self.list_of_unaligned_fasta_paths = []
        self._populate_list_of_unaligned_fasta_paths()
        self.guidance_perl_script_full_path = None
        self._get_full_path_to_guidance_perl_script()

    def _get_full_path_to_guidance_perl_script(self):
        proc = subprocess.run(['which', 'guidance'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        guidance_path = proc.stdout.decode("utf-8").rstrip()
        with open(guidance_path, 'r') as f:
            for line in f:
                if 'guidance.pl' in line:
                    self.guidance_perl_script_full_path = line.split(' ')[1]


    def _populate_list_of_unaligned_fasta_paths(self):
        # browse the output dir and get a list of all of the fastas that need to be aligned
        list_of_orth_group_dirs = list(os.walk(self.base_alignment_dir))[0][1]
        print("Populating the mp queue")
        for group_dir_name in list_of_orth_group_dirs:
            sys.stdout.write(f"\r{group_dir_name}")
            full_path = os.path.join(self.base_alignment_dir, group_dir_name)
            files_in_group_path = list(os.walk(full_path))[0][2]
            count = 0
            unaligned_file_name = None
            for file in files_in_group_path:
                if "unaligned" in file:
                    count += 1
                    unaligned_file_name = file
            if count == 0:
                raise RuntimeError(f"No unaligned file found for: {group_dir_name}")
            if count > 1:
                raise RuntimeError(f"Multiple file containing unaligned found in : {group_dir_name}")
            unaligned_file_path = os.path.join(full_path, unaligned_file_name)
            self.mp_queue.put(unaligned_file_path)


    def run_guidance_worker(self):

        # put in the STOPs
        for _ in range(self.threads):
            self.mp_queue.put('STOP')

        # a queue for putting the orth groups that failed the guidance analysis to be retried
        self.retry_queue = Queue()

        all_processes = []

        # Then start the workers
        for _ in range(self.threads):
            p = Process(target=self._guidance_worker, args=(self.mp_queue,self.retry_queue))
            all_processes.append(p)
            p.start()

        for p in all_processes:
            p.join()
    
    def make_summary_file(self):
        with open(os.path.join(self.base_alignment_dir, "guidance_aligned_fastas_summary.readme"), 'w') as f:
            f.write('DONE\n')
    
    def _guidance_worker(self, input_queue, output_queue):
        for input_fasta_path in iter(input_queue.get, 'STOP'):
            sys.stdout.write(f'\rPerforming Guidance alignment for {ntpath.basename(input_fasta_path)}')
            gw = GuidanceWorker(input_fasta_path, self.guidance_perl_script_full_path)
            if not os.path.exists(gw.aa_aligned_and_cropped_path):
                try:
                    gw._exe_guidance_analysis()
                    gw._drop_low_qual_cols_and_crop()
                except:
                    # delete the guidance output folder and add this to a list to retry
                    shutil.rmtree(gw.output_dir)
                    output_queue.put(gw.orth_group_id)
            else:
                print(f'\n{gw.aa_aligned_and_cropped_path} already exists, skipping this orth_group')


class GuidanceWorker:
    def __init__(self, input_fasta_path, guidance_perl_script_path):
        self.input_fasta_path = input_fasta_path
        self.seq_names = None
        self._get_seq_names_from_fasta()
        self.orth_group_id = ntpath.basename(self.input_fasta_path).split('_')[0]
        self.output_dir = os.path.join(os.path.dirname(self.input_fasta_path), '_'.join([self.orth_group_id, 'guidance_output']))
        self.guidance_perl_script_path = guidance_perl_script_path
        self.aa_cols_score_file_path = os.path.join(self.output_dir, f'{self.orth_group_id}.MAFFT.Guidance2_res_pair_res.PROT.scr')
        self.aa_alignment_file_path = os.path.join(self.output_dir, f'{self.orth_group_id}.MAFFT.PROT.aln')
        self.cds_alignment_file_path = os.path.join(self.output_dir, f'{self.orth_group_id}.MAFFT.aln')
        # Dataframes that will hold the aa and cds alignments that will be used for cropping and col dropping
        self.aa_df = None
        self.cds_df = None
        self.aa_cols_to_remove = []
        self.cds_cols_to_remove = []
        # output paths for the aligned and cropped aa and cds alignments
        self.aa_aligned_and_cropped_path = os.path.join(
            os.path.dirname(self.output_dir), f'{self.orth_group_id}.cropped_aligned_aa.fasta')
        self.cds_aligned_and_cropped_path = os.path.join(
            os.path.dirname(self.output_dir), f'{self.orth_group_id}.cropped_aligned_cds.fasta')
        
    def _exe_guidance_analysis(self):
        # run the guidance analysis
        # subprocess.run(
        #     ['guidance', '--seqFile', self.input_fasta_path, '--msaProgram', 'MAFFT',
        #      '--seqType', 'codon', '--outDir', self.output_fasta_path, '--bootstraps', '10',
        #      '--outOrder', 'as_input', '--colCutoff', '0.6', '--dataset', self.orth_group_id]
        #     , stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # if os.path.exists(self.output_dir):
        #     print(f"{self.output_dir} already exists, skipping guidance analysis")
        #     return
        print(['perl', self.guidance_perl_script_path, '--seqFile', self.input_fasta_path, '--msaProgram', 'MAFFT',
             '--seqType', 'codon', '--outDir', self.output_dir, '--bootstraps', '10',
             '--outOrder', 'as_input', '--colCutoff', '0.6', '--dataset', self.orth_group_id, '--proc_num', '1'])
        subprocess.run(
            ['perl', self.guidance_perl_script_path, '--seqFile', self.input_fasta_path, '--msaProgram', 'MAFFT',
             '--seqType', 'codon', '--outDir', self.output_dir, '--bootstraps', '10',
             '--outOrder', 'as_input', '--colCutoff', '0.6', '--dataset', self.orth_group_id, '--proc_num', '1'])
        foo = 'bar'

    def _get_seq_names_from_fasta(self):
        with open(self.input_fasta_path, 'r') as f:
            fasta_file = [line.rstrip() for line in f]
        self.seq_names = [fasta_file[i][1:] for i in range(0,8,2)]

    def _drop_low_qual_cols_and_crop(self):
        """At this point we have performed the guidance analysis.
        We now read in the outputs that tell us the scores of the aa and cds alignment columns
        we then use these to drop columns that have a low score.
        We also then crop the alignment in the begining and the end so that
        any columns for which one of th sequences did not have a sequence are dropped."""

        self._get_aa_cols_to_remove()

        # here we have a 0 based list of the columns that need to be dropped from the aa alignment
        # convert this to a 0 based index of the columns that need to be dropped from the cds alignment
        self._get_cds_cols_to_remove()

        # here we have the indices that need to be dropped for the cds and aa alignments
        # now read in the alignments as 2D lists, convert to pandas dataframe and perform the columns drops

        # aa first
        self._make_aa_df()

        # now drop the columns
        self._drop_cols_from_aa_df()

        # cds second
        self._make_cds_df()

        # now drop the columns
        self._drop_cols_from_cds_df()

        # here we have the cds and aa dfs that we can now do the cropping with and then finally write back out as
        # fasta files
        # go from either end deleting any columns that have a gap

        # aa first
        self.aa_df = self._crop_fasta_df(aligned_fasta_as_pandas_df_to_crop=self.aa_df)

        # cds second
        self.cds_df = self._crop_fasta_df(aligned_fasta_as_pandas_df_to_crop=self.cds_df)

        # now we just need to write out dfs
        self._write_out_aa_df()

        self._write_out_cds_df()

        # here we have the cds and the aa alignments cropped and written out
        # we can then use these as input into CODEML and the BUSTED programs

    def _write_out_cds_df(self):
        cds_fasta_list = self.pandas_df_to_fasta(pd_df=self.cds_df, seq_names=self.seq_names)
        with open(self.cds_aligned_and_cropped_path, 'w') as f:
            for line in cds_fasta_list:
                f.write(f'{line}\n')

    def _write_out_aa_df(self):
        aa_fasta_list = self.pandas_df_to_fasta(pd_df=self.aa_df, seq_names=self.seq_names)
        with open(self.aa_aligned_and_cropped_path, 'w') as f:
            for line in aa_fasta_list:
                f.write(f'{line}\n')

    def _drop_cols_from_cds_df(self):
        columns_to_keep_cds = [col for col in list(self.cds_df) if col not in self.cds_cols_to_remove]
        self.cds_df = self.cds_df[columns_to_keep_cds]

    def _make_cds_df(self):
        with open(self.cds_alignment_file_path, 'r') as f:
            cds_alignment_file_list = self.convert_interleaved_to_sequencial_fasta_two([line.rstrip() for line in f])
        self.cds_df = pd.DataFrame([list(line) for line in cds_alignment_file_list if not line.startswith('>')])

    def _drop_cols_from_aa_df(self):
        columns_to_keep_aa = [col for col in list(self.aa_df) if col not in self.aa_cols_to_remove]
        self.aa_df = self.aa_df[columns_to_keep_aa]

    def _make_aa_df(self):
        with open(self.aa_alignment_file_path, 'r') as f:
            aa_alignment_file_list = self.convert_interleaved_to_sequencial_fasta_two([line.rstrip() for line in f])
        self.aa_df = pd.DataFrame([list(line) for line in aa_alignment_file_list if not line.startswith('>')])

    def _get_cds_cols_to_remove(self):
        for col_index in self.aa_cols_to_remove:
            self.cds_cols_to_remove.extend([col_index * 3 + n for n in range(3)])

    def _get_aa_cols_to_remove(self):
        with open(self.aa_cols_score_file_path, 'r') as f:
            aa_col_score_file_list = [line.rstrip() for line in f][1:]
        for i in range(int(len(aa_col_score_file_list) / 4)):
            # the range that represents i and the next three i s in the iteration without increasing i
            file_line_indices = [i * 4 + k for k in range(4)]
            scores_set = [
                float(aa_col_score_file_list[j].split()[2]) if '-nan' not in aa_col_score_file_list[j] else '-nan' for j
                in file_line_indices]

            # now examine what is in the score sets
            drop = False
            nan_count = 0
            for score in scores_set:
                if score != '-nan':
                    if score < 0.6:
                        drop = True
                        break
                elif score == '-nan':
                    nan_count += 1
                else:
                    continue
            # when we come out ned to check if the nan_score == 4
            if nan_count == 4:
                drop = True

            # if drop = True then this is a column to drop
            if drop:
                self.aa_cols_to_remove.append(i)

    def convert_interleaved_to_sequencial_fasta_two(self, fasta_in):
        fasta_out = []

        for i in range(len(fasta_in)):

            if fasta_in[i].startswith('>'):
                if fasta_out:
                    # if the fasta is not empty then this is not the first
                    fasta_out.append(temp_seq_str)
                # else then this is the first sequence and there is no need to add the seq.
                temp_seq_str = ''
                fasta_out.append(fasta_in[i])
            else:
                temp_seq_str = temp_seq_str + fasta_in[i]
        # finally we need to add in the last sequence
        fasta_out.append(temp_seq_str)
        return fasta_out

    def _crop_fasta_df(self, aligned_fasta_as_pandas_df_to_crop):
        columns_to_drop = []
        for i in list(aligned_fasta_as_pandas_df_to_crop):
            # if there is a gap in the column at the beginning
            if '-' in list(aligned_fasta_as_pandas_df_to_crop[i]) or '*' in list(aligned_fasta_as_pandas_df_to_crop[i]):
                columns_to_drop.append(i)
            else:
                break
        for i in reversed(list(aligned_fasta_as_pandas_df_to_crop)):
            # if there is a gap in the column at the end
            if '-' in list(aligned_fasta_as_pandas_df_to_crop[i]) or '*' in list(aligned_fasta_as_pandas_df_to_crop[i]):
                columns_to_drop.append(i)
            else:
                break

        # get a list that is the columns indices that we want to keep
        col_to_keep = [col_index for col_index in list(aligned_fasta_as_pandas_df_to_crop) if
                       col_index not in columns_to_drop]
        # drop the gap columns
        return aligned_fasta_as_pandas_df_to_crop[col_to_keep]

    def pandas_df_to_fasta(self, pd_df, seq_names):
        temp_fasta = []
        for i in range(len(seq_names)):
            temp_fasta.extend(['>{}'.format(seq_names[i]), ''.join(list(pd_df.loc[i]))])
        return temp_fasta

ga = GAlign()
ga.run_guidance_worker()
ga.make_summary_file()
