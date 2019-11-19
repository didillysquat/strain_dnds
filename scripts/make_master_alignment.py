from collections import defaultdict
import sys
import os
import subprocess
class MasterAlignment:
    def __init__(self):
        self.species = sys.argv[1]
        self.base_input_dir = os.path.join(
            '/home/humebc/projects/parky/breviolum_transcriptomes/local_alignments', self.species)
        self.list_of_prottest_paths = []
        self.list_of_prottest_files_not_found = []
        self._populate_list_of_prottest_paths()
        self.master_fasta_out_path = sys.argv[2]
        os.makedirs(os.path.dirname(self.master_fasta_out_path), exist_ok=True)
        self.q_file_out_path = sys.argv[3]
        os.makedirs(os.path.dirname(self.q_file_out_path), exist_ok=True)
        self.bad_prot_file_list = []
        self.prottest_path = '/home/humebc/phylogeneticSoftware/prottest/prottest-3.4.2/prottest-3.4.2.jar'
        
    def _populate_list_of_prottest_paths(self):
        list_of_orth_group_dirs = list(os.walk(self.base_input_dir))[0][1]
        print("Populating the list of prottest paths")
        for group_dir_name in list_of_orth_group_dirs:
            sys.stdout.write(f"\r{group_dir_name}")
            path_to_prot_out_file = os.path.join(self.base_input_dir, group_dir_name, f'{group_dir_name}_prottest_result.out')
            if os.path.exists(path_to_prot_out_file):
                self.list_of_prottest_paths.append(path_to_prot_out_file)
            else:
                self.list_of_prottest_files_not_found.append(path_to_prot_out_file)
        
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
            for j in range(300, len(temp_file_list), 1):
                if 'Best model according to BIC' in temp_file_list[j]:
                    model = temp_file_list[j].split(':')[1].strip().replace('+G', '').replace('+I', '')
                    break
            if model == '':
                # attempt to rerun the protein model once
                os.remove(protpath)
                subprocess.run(
                    ['java', '-jar', self.prottest_path, '-i',
                     protpath.replace('_prottest_result.out', '.cropped_aligned_aa.fasta'), '-o', protpath,
                     '-all-distributions', '-all'])
                model = ''
                # orth_num = int(os.path.dirname(protpath).split('_')[0])
                with open(protpath, 'r') as f:
                    temp_file_list = [line.rstrip() for line in f]
                for j in range(300, len(temp_file_list), 1):
                    if 'Best model according to BIC' in temp_file_list[j]:
                        model = temp_file_list[j].split(':')[1].strip().replace('+G', '').replace('+I', '')
                        break
                if model == '':
                    self.bad_prot_file_list.append(protpath)
                    continue
                # sys.exit('Model line not found in {}'.format(protpath))
            model_to_orth_dict[model].append(protpath.replace('_prottest_result.out', '.cropped_aligned_aa.fasta'))

        # #N.B. that we cannot have different gamma for different partitions
        # # Also best advice is not to run +G and +I together.
        # # As such we only need to extract the base model here i.e. WAG rather than WAG [+G|+I]
        # for model in model_to_orth_dict

        # delete the bad prot models and rerun the prottest to see if we can get them good.
        # for bad_path in self.bad_prot_file_list:
        #     os.remove(bad_path)

        print('The 19k sequences are best represented by {} different aa models'.format(len(model_to_orth_dict.keys())))

        # here we have the dict populated
        # now sort the dict
        sorted_model_list = sorted(model_to_orth_dict, key=lambda k: len(model_to_orth_dict[k]), reverse=True)

        # now go model by model in the sorted_model_list to make the master alignment.

        # not the most elegant way but I think I'll just create the master fasta in memory
        if self.species == 'b_minutum':
            master_fasta = ['>SRR1793320', '', '>SRR1793321', '', '>SRR1793322', '', '>SRR1793323', '']
        elif self.species == 'b_psygmophilum':
            master_fasta = ['>SRR1793324', '', '>SRR1793325', '', '>SRR1793326', '', '>SRR1793327', '']

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
        with open(self.master_fasta_out_path, 'w') as f:
            for line in master_fasta:
                f.write('{}\n'.format(line))

        # now write out the q file
        with open(self.q_file_out_path, 'w') as f:
            for line in q_file:
                f.write('{}\n'.format(line))

ma = MasterAlignment()
ma.concatenate_local_alignment()