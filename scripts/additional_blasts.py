"""This script will read in a blast output from a single blast of a set of transcripts against
the swiss prot database. It will also read in the long isoform only version of the assembly fasta.
If there are transcripts that have not found a match in this blast output base, they will in turn be blasted against
the trembl database. Any non-matches after this blast will then be blast against the ncbi_nr database.
The ouputs of this will be a blast match for each of the trembl and ncbi_nr database.

paths are passed to this script in the following order:
{input.swiss_prot_blast_output} {input.trembl_db_path_string} {input.long_iso_only_transcript_fasta}
{input.ncbi_db_path_string} {output.trembl_blast_output} {output.ncbi_nr_blast_output}
{output.consolidated_blast_out_results}

arg[1] = the path to the blast match file for the swiss prot db blast
arg[2] = the path to the long isoform only transcript fasta
arg[3] = the string that should be provided to the -db argument when running the trembl blast
arg[4] = the string that should be provided to the -db argument when running the ncbi_nr blast
arg[5] = the output file path for the trembl blast (provided to -out for the blast)
arg[6] = the output file path for the ncbi_nr blast (provided to -out for the blast)
arg[7] = the file path for the consolidated blast output results
arg[8] = the file path for the tremble fasta of no match sequences to be input to the trembl blast
arg[9] = the file path for the ncbi_nr fasta of no match sequences to be input to the ncbi_nr blast
arg[10] = number of threads
arg[11] = species
arg[12] = sra
"""

from Bio import SeqIO
import sys
from collections import defaultdict
import subprocess
import os

class additionalBlast():
    def __init__(self):
        # get the trinity_fasta
        self.trinity_fasta = list(SeqIO.parse(sys.argv[2], "fasta"))
        # The transcript IDs that have matches
        self.matched_seqs_list = set()
        # A new trinity_fasta with only the unmatched seqs
        self.unmatched_fasta = []
        # A new version of the blast results with only one match per query sequence (best match kept)
        # key will be the query seq name, value will be the full string of the match as in the blast output
        self.best_match_collection = {}
        # swissprot output blast file path
        self.swiss_prot_blast_out_file_path = sys.argv[1]
        # tremble output blast file path
        self.trembl_blast_out_file_path = sys.argv[5]
        # ncbi_pr output blast file path
        self.ncbi_nr_blast_out_file_path = sys.argv[6]
        # the output path for the consolidated blast output
        self.consolidated_blast_output_path = sys.argv[7]
        # the path to use for the trembl input fasta
        self.tremble_fasta_in_path = sys.argv[8]
        # the path to use for the ncbi_nr input fasta
        self.ncbi_nr_fasta_in_path = sys.argv[9]
        # the string path to the trembl db
        self.trembl_db_path = sys.argv[3]
        # the string path to the ncbi_nr db
        self.ncbi_db_path = sys.argv[4]
        # number of threads to use for the blasts
        self.threads = sys.argv[10]
        self.species = sys.argv[11]
        self.sra = sys.argv[12]
        # check that all directories exist and make them if not
        os.makedirs(os.path.dirname(self.trembl_blast_out_file_path), exist_ok=True)
        os.makedirs(os.path.dirname(self.ncbi_nr_blast_out_file_path), exist_ok=True)

    def do_additional_blasts(self):
        self.check_blast(blast_id='swiss_prot', output_file=self.swiss_prot_blast_out_file_path)

    def check_blast(self, output_file, blast_id='swiss_prot'):
        # go through the matches and look to see which sequence have a match and create the best_match_collection,
        # as well as the matched_seqs_list and the self.unmatched_fasta. If the unmatched fasta has items in it
        # then we will need to write out the fasta and perform a blast.
        # Else simply write out empty blast results and make the consolidated output

        print(f'Checking the {blast_id} blast output results for {self.species} {self.sra}')
        print(f"Curating best matches for {blast_id}_{self.species}_{self.sra}")
        # get the swiss prot blast matches
        with open(output_file, 'r') as f:
            blast_matches = f.readlines()
        # A dict to keep track of the best match for a given query sequence
        best_match_dd = defaultdict(int)
        for blast_match_line in blast_matches:
            split_list = blast_match_line.split('\t')
            seq_id = split_list[0]
            bit_score = float(split_list[11])
            if best_match_dd[seq_id]:
                if best_match_dd[seq_id] < bit_score: # current match is a better match
                    best_match_dd[seq_id] = bit_score
                    self.best_match_collection = blast_match_line

                else:
                    pass # This is not as good a match as a previous match
            else:
                # Then there is not yet a best match for this sequence
                best_match_dd[seq_id] = bit_score
                self.best_match_collection = blast_match_line
                self.matched_seqs_list.add(seq_id)


        # Here we have the best_match_collection populated and we can now compare this to the trinity_fasta
        # to see which of the sequence ids still remain unmatched.
        # An example seq_id from the blast match file: TRINITY_DN37624_c0_g1_i1
        self.unmatched_fasta = [] # Important to reset this between blast ouput checks
        for seq_rec in self.trinity_fasta:
            if not seq_rec.name in self.matched_seqs_list:
                self.unmatched_fasta.append(seq_rec)

        print(f"{len(self.unmatched_fasta)} unmatched seq_ids remaining after the {blast_id} blast for {self.species} {self.sra}")
        # here we have a list of seq records for thos transcripts that have not been matched using the
        # swiss prot database.
        # These should be blasted against the trembl database
        if self.unmatched_fasta:
            if blast_id == 'swiss_prot':
                print(f"Performing tremble blast for {self.species} {self.sra}")
                self.do_blast(fasta_in_path=self.tremble_fasta_in_path, blast_out_path=self.trembl_blast_out_file_path, db_path=self.trembl_db_path)
                self.check_blast(blast_id='trembl', output_file=self.trembl_blast_out_file_path)
            elif blast_id == 'trembl':
                print(f"Performing ncbi_nr blast for {self.species} {self.sra}")
                self.do_blast(fasta_in_path=self.ncbi_nr_fasta_in_path, blast_out_path=self.ncbi_nr_blast_out_file_path, db_path=self.ncbi_db_path)
                self.check_blast(blast_id='ncbi_nr', output_file=self.ncbi_nr_blast_out_file_path)
        print(f"Outputting consolidated blast results for {self.species} {self.sra}")
        self.finalise_results_and_exit()

    def do_blast(self, fasta_in_path, blast_out_path, db_path):
        # write out the no match fasta
        # write it out manual so that the fasta is not interleaved
        with open(fasta_in_path, 'w') as f:
            for seq_rec in self.unmatched_fasta:
                f.write(f'>{seq_rec.description}\n')
                f.write(f'{seq_rec._seq}\n')

        # now do the blast
        subprocess.run(
            ['blastx', '-out', blast_out_path, '-outfmt', '6', '-query', fasta_in_path, '-db', db_path,
             '-evalue', '1e-5', '-num_threads', self.threads])

    def finalise_results_and_exit(self):
        """Write out the consolidated blast result file using the self.best_match_collection."""
        # now write out the best_match_collection
        with open(self.consolidated_blast_output_path, 'w') as f:
            for blast_line in self.best_match_collection.values():
                f.write(blast_line)
        return

ab = additionalBlast()
ab.do_additional_blasts()