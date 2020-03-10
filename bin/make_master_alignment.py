#!/usr/bin/env python3
from collections import defaultdict
import sys
import os
import subprocess
class MasterAlignment:
    def __init__(self):
        self.prot_out_dir = sys.argv[1]
        self.list_of_prottest_paths = self._populate_list_of_prottest_paths()
        self.bad_prot_file_list = []
        self.good_prot_out_count = 0
        
    def _populate_list_of_prottest_paths(self):
        for root, dirs, files in os.walk(self.prot_out_dir):
            return [os.path.join(self.prot_out_dir, file_name) for file_name in files if '.out' in file_name]
        
    def concatenate_local_alignment(self):
        # The master alignment that we create should be partitioned according to the protein model used.
        # I have generated all of the .out files which are the outputs from the prottest
        # We should iter through these and create a dictionary that is a model type as key
        # and then have a list of the orthologs of that model.
        # then sort this by the length of the list
        # then work our way through the local alignments in this order creating the alignment
        # We will need to generate the p file as we go
        # this should take the form
        '''
        JTT, gene1 = 1-500
        WAGF, gene2 = 501-800
        WAG, gene3 = 801-1000
        '''

        # iter through the list of protfiles creating the dict relating model to ortholog
        # we cannnot change the +G or +I for each partition. As such I will define according to the base model
        model_to_orth_dict = defaultdict(list)
        for protpath in self.list_of_prottest_paths:
            sys.stdout.write(f'\r{protpath}')
            model = ''
            # orth_num = int(os.path.dirname(protpath).split('_')[0])
            with open(protpath, 'r') as f:
                temp_file_list = [line.rstrip() for line in f]
            if not temp_file_list:
                continue
            for j in range(0, len(temp_file_list), 1):
                if 'Best model according to BIC' in temp_file_list[j]:
                    if '+' in temp_file_list[j+2]:
                        model = temp_file_list[j+2].split(':')[1].strip().split('+')[0]
                    else:
                        model = temp_file_list[j+2].split(':')[1].strip()
                    if '-' in model:
                        model = model.split('-')[0]
                    break
            if model == '':
                # # attempt to rerun the protein model once
                # os.remove(protpath)
                # subprocess.run(
                #     ['java', '-jar', self.prottest_path, '-i',
                #      protpath.replace('_prottest_result.out', '.cropped_aligned_aa.fasta'), '-o', protpath,
                #      '-all-distributions', '-all'])
                # model = ''
                # # orth_num = int(os.path.dirname(protpath).split('_')[0])
                # with open(protpath, 'r') as f:
                #     temp_file_list = [line.rstrip() for line in f]
                # for j in range(300, len(temp_file_list), 1):
                #     if 'Best model according to BIC' in temp_file_list[j]:
                #         model = temp_file_list[j].split(':')[1].strip().replace('+G', '').replace('+I', '')
                #         break
                # if model == '':
                #     self.bad_prot_file_list.append(protpath)
                #     continue
                # sys.exit('Model line not found in {}'.format(protpath))
                continue
            self.good_prot_out_count += 1
            model_to_orth_dict[model].append(protpath.replace('_prottest_result.out', '_cropped_aligned_aa.fasta'))

        # #N.B. that we cannot have different gamma for different partitions
        # # Also best advice is not to run +G and +I together.
        # # As such we only need to extract the base model here i.e. WAG rather than WAG [+G|+I]
        # for model in model_to_orth_dict

        print(f'The {self.good_prot_out_count} sequences are best represented by {len(model_to_orth_dict.keys())} different aa models')

        # here we have the dict populated
        # now sort the dict
        sorted_model_list = sorted(model_to_orth_dict, key=lambda k: len(model_to_orth_dict[k]), reverse=True)

        # now go model by model in the sorted_model_list to make the master alignment.

        # not the most elegant way but I think I'll just create the master fasta in memory
        # get a list of the strain names to use in making the master fasta
        an_aa_alignment_path = model_to_orth_dict[sorted_model_list[0]][0]
        with open(an_aa_alignment_path, 'r') as f:
            lines_list = [line.rstrip() for line in f]
        # The below should produce e.g. ['>SRR1793324', '', '>SRR1793325', '', '>SRR1793326', '', '>SRR1793327', ''] (for four strains)
        master_fasta = [f">{'_'.join(lines_list[i].split('_')[1:])}" if (i%2 == 0) else '' for i in range(len(lines_list))]
        # master_fasta = ['>BreviolumB5_Sradians', '', '>Breviolumfaviinorum_Pclivosa', '', '>SRR1793324', '', '>SRR1793325', '', '>SRR1793326', '', '>SRR1793327', '']

        # The q file will hold the information for the partitioning of the alignment for the raxml analysis
        q_file = []
        for model in sorted_model_list:
            q_file_start = len(master_fasta[1]) + 1
            sys.stdout.write('\rProcessing model {} sequences'.format(model))
            for aa_align_path in model_to_orth_dict[model]:
                with open(aa_align_path, 'r') as f:
                    temp_list_of_lines = [line.rstrip() for line in f]

                for i in range(1, len(temp_list_of_lines), 2):
                    new_seq = master_fasta[i] + temp_list_of_lines[i]
                    master_fasta[i] = new_seq
            q_file_finish = len(master_fasta[1])
            q_file.append('{}, gene{} = {}-{}'.format(
                model.upper(), sorted_model_list.index(model) + 1, q_file_start, q_file_finish))

        # here we have the master fasta and the q file ready to be outputted

        # now write out the master fasta
        with open("master_fasta_for_tree.fasta", 'w') as f:
            for line in master_fasta:
                f.write('{}\n'.format(line))

        # now write out the q file
        with open("q_partition_file.q", 'w') as f:
            for line in q_file:
                f.write('{}\n'.format(line))

ma = MasterAlignment()
ma.concatenate_local_alignment()