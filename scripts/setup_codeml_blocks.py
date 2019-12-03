import os
import sys
import pickle

class CODEML:
    def __init__(self):
        self.species = sys.argv[1]
        self.base_input_dir = os.path.join(
            '/home/humebc/projects/parky/breviolum_transcriptomes/local_alignments', self.species)
        self.output_dir = os.path.join('/home/humebc/projects/parky/breviolum_transcriptomes/codeml', self.species)
        os.makedirs(self.output_dir, exist_ok=True)
        self.list_of_orth_group_dirs = list(os.walk(self.base_input_dir))[0][1]
        self.bad_dir_list = []
        # counter that we'll use to split into 20 roughly 1k sequence alignments
        self.counter = 0
        # counter for keeping track of which block we are on
        self.block_counter = 0
        # This is a list of the tuples that contain the paths to the blocks and info needed
        # to run the codeml analysis
        self.list_of_guidance_dirs = []
        self.block_info_pickle_out_path = sys.argv[2]
        self.phylip_alignment = []


    def generate_block_phylip_alignments_for_CODEML(self):
        ''' The documentation for what format the files should be in for submitting to PAML/CODEML are not so
        great. But from some testing and looking through the examples that are available two things seem to be key
        1 - automatic pairwise comparisons of sequences can be performed using the runmode as -2
        2 - a blocked alignment format is the easiest way to test all of the orthologs at once but unlike
        in the documentation, (using the G marker) this is easiest done by simply having a phymil alignment format
        one after another in succestion in a sinlge file and then setting the ndata setting to how ever many
        sets of aligments you have (about 19K for us). I've run some tests and this all seems to be running well.
        Unfortunately, the only way to MP this seems to be physically splitting up the data and running multiple
        instance of the CODEML executable.
        So... with this in mind, I will split up the 19K or so sequences into phylip alignments of 1000 sequences.
        I will then make respective control files for each of these, put them in their own directory and then set
        an instance of the CODEML running for each of these.'''

        # for each directory or ortholog
        tot_num_dirs = len(self.list_of_orth_group_dirs)
        from_bad = False
        for orth_dir in self.list_of_orth_group_dirs:
            sys.stdout.write(f'\rOrtholog: {orth_dir}     {self.counter}/{tot_num_dirs}')

            # make a new alignment file for every 1000 individual alignments
            if self.counter % 1000 == 0 and not from_bad:
                if self.block_counter != 0:
                    # then we already have a block that needs writing
                    seq_file_ctrl_file_tup = self.write_out_cntrl_and_seq_file(
                        phylip_alignment=self.phylip_alignment, num_align=1000)
                    self.list_of_guidance_dirs.append(seq_file_ctrl_file_tup)
                # once the old block is written out start a new one
                self.block_counter += 1
                if self.block_counter == 16:
                    foo = 'asdf'
                os.makedirs(f'{self.output_dir}/block_{self.block_counter}', exist_ok=True)
                print(f'\n\nStarting block {self.block_counter}')
                self.phylip_alignment = []

            # add single phylip alignment to block master alignment
            # if the fasta was empty then this will return False
            single_phylip = self.generate_phylip_from_fasta(orth_dir)

            if single_phylip:
                self.phylip_alignment.extend(single_phylip)
                self.counter += 1
                from_bad = False
            else:
                # if the fasta was empty then log this and don't add anything to the counter
                self.bad_dir_list.append(orth_dir)
                from_bad = True

        # now write out the final block of alignments
        seq_file_ctrl_file_tup = self.write_out_cntrl_and_seq_file(
            phylip_alignment=self.phylip_alignment,
            num_align=len(self.list_of_orth_group_dirs)-len(self.bad_dir_list)-(1000*(self.block_counter-1)))
        self.list_of_guidance_dirs.append(seq_file_ctrl_file_tup)
        pickle.dump(self.list_of_guidance_dirs, open(self.block_info_pickle_out_path, "wb" ))

    def generate_phylip_from_fasta(self, orth_grp_id):
        temp_str = str()
        temp_list = list()
        cds_fasta_path = os.path.join(self.base_input_dir, orth_grp_id, f'{orth_grp_id}.cropped_aligned_cds.fasta')
        if not os.path.exists(cds_fasta_path):
            return False

        with open(cds_fasta_path, 'r') as f:
            fasta_file = [line.rstrip() for line in f]

        if len(fasta_file[1]) == 0:
            # if the fasta is empty then we need to log this outside
            return False
        else:
            for line in fasta_file:
                if line.startswith('>'):
                    temp_str = line[1:]
                else:
                    temp_list.append('{}    {}'.format(temp_str, line))
            # finally put in the header file
            temp_list.insert(0, f'\t{len(temp_list)} {len(fasta_file[1])}')

            return temp_list


    def write_out_cntrl_and_seq_file(self, phylip_alignment, num_align):
        # write out the control file specific to this alignment
        ctrl_file_path = self.write_out_control_file(num_alignments=num_align)
        # write out the phylip file
        seq_file_path = os.path.join(
            self.output_dir, f'block_{self.block_counter}', f'block_{self.block_counter}_cds.phylip')
        with open(seq_file_path, 'w') as f:
            for line in phylip_alignment:
                f.write('{}\n'.format(line))
        # write out the tree file
        tree_file = None
        if self.species == 'b_minutum':
            tree_file = '((SRR1793323:0.00207691863248198353,SRR1793321:0.00135412970390987497):' \
                        '0.00162401761552979072,SRR1793322:0.00313910729521391747,' \
                        'SRR1793320:0.00302742649625364728):0.0;'
        elif self.species == 'b_psygmophilum':
            tree_file = '(SRR1793326:0.00128083627808537638,(SRR1793327:0.00133413967614412271,SRR1793325:' \
                        '0.00122248010329080036):0.00046190387263074311,SRR1793324:0.00327381204178443319):0.0;'

        tree_file_path = os.path.join(self.output_dir, f'block_{self.block_counter}', f'block_{self.block_counter}_tree.nwk')
        with open(tree_file_path, 'w') as f:
            f.write(f'{tree_file}\n')

        return (seq_file_path, ctrl_file_path, tree_file_path)

    def write_out_control_file(self, num_alignments):
        block_output_dir = os.path.join(self.output_dir, f'block_{self.block_counter}')
        seq_file_path = os.path.join(block_output_dir, f'block_{self.block_counter}_cds.phylip')
        out_file_path = os.path.join(block_output_dir, f'block_{self.block_counter}_codeml_results.out')
        ctrl_path     = os.path.join(block_output_dir, f'block_{self.block_counter}_cds.ctrl')
        tree_file_path = os.path.join(block_output_dir, f'block_{self.block_counter}_tree.nwk')
        ctrl_file = [
        f'seqfile = {seq_file_path}',
        f'treefile = {tree_file_path}',
        f'outfile = {out_file_path}',
        'runmode = -2',
        'seqtype = 1',
        'CodonFreq = 2',
        f'ndata = {num_alignments}',
        'clock = 0',
        'model = 0',
        'NSsites = 0',
        'icode = 0',
        'fix_omega = 0',
        'omega = .4',
        'cleandata = 0',
        'verbose = 1'
        ]

        with open(ctrl_path, 'w') as f:
            for line in ctrl_file:
                f.write('{}\n'.format(line))

        return ctrl_path

cml = CODEML()
cml.generate_block_phylip_alignments_for_CODEML()