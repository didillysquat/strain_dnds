import os
import sys
import pickle

class CODEML:
    def __init__(self):
        self.species = sys.argv[1]
        self.base_input_dir = os.path.join(
            '/home/humebc/projects/parky/breviolum_transcriptomes/local_alignments', self.species)
        self.output_dir = '/home/humebc/projects/parky/breviolum_transcriptomes/local_alignments'
        os.makedirs(self.output_dir, exist_ok=True)
        self.list_of_orth_group_dirs = list(os.walk(self.base_input_dir))[0][1]
        self.bad_dir_list = []
        # counter that we'll use to split into 20 roughly 1k sequence alignments
        self.counter = 0
        # counter for keeping track of which block we are on
        self.block_counter = 0
        # I dont know what this is used for yet
        self.list_of_guidance_dirs = []
        self.summary_file_path = sys.argv[1]
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
        for orth_dir in self.list_of_orth_group_dirs:
            sys.stdout.write(f'\rOrtholog: {orth_dir}     {self.counter}/{tot_num_dirs}')

            # make a new alignment file for every 1000 individual alignments
            if self.counter % 1000 == 0:
                if self.block_counter != 0:
                    # then we already have a block that needs writing
                    seq_file_ctrl_file_tup = self.write_out_cntrl_and_seq_file(
                        block_counter=self.block_counter, output_dir=self.output_dir,
                        phylip_alignment=self.phylip_alignment, num_align=1000)
                    self.list_of_guidance_dirs.append(seq_file_ctrl_file_tup)
                # once the old block is written out start a new one
                self.block_counter += 1
                os.makedirs(f'{self.output_dir}/block_{self.block_counter}', exist_ok=True)
                print(f'\n\nStarting block {self.block_counter}')
                self.phylip_alignment = []

            # add single phylip alignment to block master alignment
            # if the fasta was empty then this will return False
            single_phylip = self.generate_phylip_from_fasta(orth_dir)

            if single_phylip:
                self.phylip_alignment.extend(single_phylip)
                self.counter += 1
            else:
                # if the fasta was empty then log this and don't add anything to the counter
                self.bad_dir_list.append(orth_dir)

        # write out a list of the poor alignement orfs
        with open(self.summary_file_path, 'w') as f:
            for line in self.bad_dir_list:
                f.write(f'{line}\n')

        # now write out the final block of alignments
        seq_file_ctrl_file_tup = self.write_out_cntrl_and_seq_file(self.block_counter, self.output_dir, self.phylip_alignment, num_align=len(self.list_of_orth_group_dirs)-len(self.bad_dir_list)-(1000*(self.block_counter-1)))
        self.list_of_guidance_dirs.append(seq_file_ctrl_file_tup)

        pickle.dump(self.list_of_guidance_dirs, open('{}/list_of_guidance_dirs.pickle'.format(self.output_dir), "wb" ))

    def generate_phylip_from_fasta(self, orth_grp_id):
        temp_str = str()
        temp_list = list()

        with open(os.path.join(self.base_input_dir, orth_grp_id, f'{orth_grp_id}.cropped_aligned_aa.fasta'), 'r') as f:
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


    def write_out_cntrl_and_seq_file(self, block_counter, output_dir, phylip_alignment, num_align):
        # write out the control file specific to this alignment
        ctrl_file_path = self.write_out_control_file(
            output_dir='{}/block_{}'.format(output_dir, block_counter),
            num_alignments=num_align,
            grp=block_counter)
        # write out the phylip file
        seq_file_path = '{}/block_{}/block_{}_cds.phylip'.format(output_dir, block_counter, block_counter)
        with open(seq_file_path, 'w') as f:
            for line in phylip_alignment:
                f.write('{}\n'.format(line))
        # write out the tree file
        tree_file = '(ppsyg:0.01524804457090833502,(min:0.00305561548082329418,pmin:0.00350296114601793013)' \
                    ':0.03350192310501232812,psyg:0.01618135662493049715);'
        tree_file_path = '{}/block_{}/block_{}_tree.nwk'.format(output_dir, block_counter, block_counter)
        with open(tree_file_path, 'w') as f:
            f.write('{}\n'.format(tree_file))

        return (seq_file_path, ctrl_file_path, tree_file_path)

    def write_out_control_file(self, output_dir, num_alignments, grp):
        seq_file_path = '{}/block_{}_cds.phylip'.format(output_dir, grp)
        out_file_path = '{}/block_{}_guidance_results.out'.format(output_dir, grp)
        ctrl_path     = '{}/block_{}_cds.ctrl'.format(output_dir, grp)
        tree_file_path = '{}/block_{}/block_{}_tree.nwk'.format(output_dir, grp, grp)
        ctrl_file = [
        'seqfile = {}'.format(seq_file_path),
        'treefile = {}'.format(tree_file_path),
        'outfile = {}'.format(out_file_path),
        'runmode = -2',
        'seqtype = 1',
        'CodonFreq = 2',
        'ndata = {}'.format(num_alignments),
        'clock = 0',
        'model = 0',
        'NSsites = 0',
        'icode = 0',
        'fix_omega = 0',
        'omega = .4',
        'cleandata = 0'
        ]

        with open(ctrl_path, 'w') as f:
            for line in ctrl_file:
                f.write('{}\n'.format(line))

        return ctrl_path