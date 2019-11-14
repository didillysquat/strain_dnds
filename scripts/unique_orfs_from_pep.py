import sys
from Bio import SeqIO
from collections import defaultdict

class unique_orfs:
    def __init__(self):
        self.pep_file = list(SeqIO.parse(sys.argv[1], "fasta"))
        self.out_pep_path = sys.argv[2]

    def unique(self):
        # dd to hold the list of orfs that have been predicted for each transcript
        orf_dd_list = defaultdict(list)


        # populate the orf_dd dict
        tot_len = len(self.pep_file)
        count = 0
        for seq_rec in self.pep_file:
            count += 1
            sys.stdout.write(f'\r{((count/tot_len)*100):.2f} percent complete')
            elements = seq_rec.description.split(' ')
            transcript_id = elements[0].split('.')[0]
            print(f'populating dd dict with {transcript_id}')
            completeness = elements[1].split(':')[1]
            if completeness == 'complete':
                complete = True
            else:
                complete = False
            length = int(elements[2].split(':')[1])
            orf_dd_list[transcript_id].append(ORF(complete=complete, length=length, record=seq_rec))


        # Now for each of the gene ids, choose the longest complete ORF
        # If there is not a complete orf just take the longest
        final_pep_list = []
        for gene_id, orf_list in orf_dd_list.items():
            longest_complete = 0
            complete_best_match = None
            longest_non_complete = 0
            non_complete_best_match = None
            for orf in orf_list:
                if orf.complete:
                    if orf.length > longest_complete:
                        complete_best_match = orf.record
                        longest_complete = orf.length
                else:
                    if orf.length > longest_non_complete:
                        non_complete_best_match = orf.record
                        longest_non_complete = orf.length

            # if a complete exists then we will use this
            # else use the non-complete
            if complete_best_match:
                final_pep_list.append(complete_best_match)
            else:
                final_pep_list.append(non_complete_best_match)

        # write out the list of sequence records manually so that they aren't interleaved
        print(f'Writing out new unique pep file to: {self.out_pep_path}')
        with open(self.out_pep_path, 'w') as f:
            for seq_rec in final_pep_list:
                f.write(f'>{seq_rec.description}\n')
                f.write(f'{seq_rec._seq}\n')



class ORF:
    def __init__(self, complete, length, record):
        self.complete = complete
        self.length = length
        self.record = record

uo = unique_orfs()
uo.unique()