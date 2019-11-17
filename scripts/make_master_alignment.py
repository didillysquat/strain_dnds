from collections import defaultdict
import sys
import os

class MasterAlignment:
    def __init__(self):
        self.species = sys.argv[1]
        self.base_input_dir = os.path.join(
            '/home/humebc/projects/parky/breviolum_transcriptomes/local_alignments', self.species)
        self.list_of_prottest_paths = []
        self._populate_list_of_prottest_paths()
        self.list_of_prottest_files_not_found = []
        self.master_fasta_out_path = sys.argv[2]
        self.q_file_out_path = sys.argv[3]
        
    def _populate_list_of_prottest_paths(self):
        list_of_orth_group_dirs = list(os.walk(self.base_input_dir))[0][1]
        print("Populating the list of prottest paths")
        for group_dir_name in list_of_orth_group_dirs:
            sys.stdout.write(f"\r{group_dir_name}")
            path_to_prot_out_file = os.path.join(group_dir_name, '_prottest_result.out')
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
            model = ''
            # orth_num = int(os.path.dirname(protpath).split('_')[0])
            with open(protpath, 'r') as f:
                temp_file_list = [line.rstrip() for line in f]
            for j in range(300, len(temp_file_list), 1):
                if 'Best model according to BIC' in temp_file_list[j]:
                    model = temp_file_list[j].split(':')[1].strip().replace('+G', '').replace('+I', '')
                    break
            if model == '':
                # sys.exit('Model line not found in {}'.format(orth_num))
                sys.exit('Model line not found in {}'.format(protpath))
            model_to_orth_dict[model].append(protpath.replace('_prottest_result.out', '.cropped_aligned_aa.fasta'))

        # #N.B. that we cannot have different gamma for different partitions
        # # Also best advice is not to run +G and +I together.
        # # As such we only need to extract the base model here i.e. WAG rather than WAG [+G|+I]
        # for model in model_to_orth_dict

        print('The 19k sequences are best represented by {} different aa models'.format(len(model_to_orth_dict.keys())))

        # here we have the dict populated
        # now sort the dict
        sorted_model_list = sorted(model_to_orth_dict, key=lambda k: len(model_to_orth_dict[k]), reverse=True)

        # now go model by model in the sorted_model_list to make the master alignment.

        # not the most elegant way but I think I'll just create the master fasta in memory
        master_fasta = ['>min', '', '>pmin', '', '>psyg', '', '>ppsyg', '']

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

        # now run raxml
        # NB note that although we are specificing mdels for each partition, we still need to add the PROTGAMMAIWAG
        # model argument to the -m flag. This is just a weird operation of raxml and is explained but hidden in the manual
        # (search for 'If you want to do a partitioned analysis of concatenated'). Raxml will only extract the rate
        # variation information from this and will ignore the model component e.g. WAG. FYI any model could be used
        # doesn't have to be WAG.
        raxml_path = '/home/humebc/phylogeneticsSoftware/raxml/standard-RAxML/raxmlHPC-PTHREADS-AVX'
        subprocess.run([raxml_path, '-s', master_fasta_output_path, '-q', q_file_output_path,
                        '-x', '183746', '-f', 'a', '-p', '83746273', '-#', '1000', '-T', '8', '-n', 'parkinson_out',
                        '-m', 'PROTGAMMAWAG', '-w', '/home/humebc/projects/parky/aa_tree_creation'])

        print('\nConstruction of master fasta complete:\n{}\n{}'.format(master_fasta_output_path, q_file_output_path))