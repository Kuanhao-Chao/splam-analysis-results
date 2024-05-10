import os
import pprint

SAMPLE_SIZE = 1000
INPUT_DIR = '../../../train/results/train_test_dataset/'
NEG_DATASETS = [f'{INPUT_DIR}input_neg_1/800bp/test_neg_1.shuffle.fa', f'{INPUT_DIR}input_neg_random/800bp/test_neg_random.shuffle.fa']
POS_DATASETS = [f'{INPUT_DIR}input_pos_alts/800bp/test_pos_alts.shuffle.fa', f'{INPUT_DIR}input_pos_mane/800bp/test_pos_mane.shuffle.fa']
OUTPUT_DIR = '../results/splam/1/'
os.makedirs(os.path.dirname(OUTPUT_DIR), exist_ok=True)

neg_donor_motifs = {}
pos_donor_motifs = {}
neg_acceptor_motifs = {}
pos_acceptor_motifs = {}

# position of donor = 200:202       -> GT
# position of acceptor = 598:600    -> AG
# just sample 1k sequences from train and then format them as input files for all of them

def reverse_complement(sequence):
    # define the complement for each nucleotide
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # replace each nucleotide by its complement and reverse the sequence
    return ''.join(complement[nuc] for nuc in reversed(dna_sequence))

def add_motif(donor_motif, acceptor_motif, group):
    if group == 'neg':
        neg_donor_motifs[donor_motif] = neg_donor_motifs.get(donor_motif, 0) + 1
        neg_acceptor_motifs[acceptor_motif] = neg_acceptor_motifs.get(acceptor_motif, 0) + 1
    if group == 'pos':
        pos_donor_motifs[donor_motif] = pos_donor_motifs.get(donor_motif, 0) + 1
        pos_acceptor_motifs[acceptor_motif] = pos_acceptor_motifs.get(acceptor_motif, 0) + 1

def display_stats():
    
    print('NEGATIVE DATASET\n', '-'*120)
    print('Donor:', sum(neg_donor_motifs.values()))
    pprint.pprint(neg_donor_motifs)
    print('Acceptor:', sum(neg_acceptor_motifs.values()))
    pprint.pprint(neg_acceptor_motifs)
    print('\nPOSITIVE DATASET\n', '-'*120)
    print('Donor:', sum(pos_donor_motifs.values()))
    pprint.pprint(pos_donor_motifs)
    print('Acceptor:', sum(pos_acceptor_motifs.values()))
    pprint.pprint(pos_acceptor_motifs)

def sample(neg_datasets, pos_datasets, neg_file, pos_file, stats_file=''):
    overall_neg = []
    overall_pos = []
    
    for dataset_file in neg_datasets:
        with open(dataset_file, 'r') as fr:
            neg_lines = [next(fr) for _ in range(2 * SAMPLE_SIZE)]
            for i, line in enumerate(neg_lines): 
                if i % 2 == 0:
                    continue
                donor_motif = line[200:202]
                acceptor_motif = line[598:600]
                add_motif(donor_motif, acceptor_motif, 'neg')

            overall_neg.extend(neg_lines)
    
    with open(neg_file, 'w') as fw: 
        fw.writelines(overall_neg)

    for dataset_file in pos_datasets:
        with open(dataset_file, 'r') as fr:
            pos_lines = [next(fr) for _ in range(2 * SAMPLE_SIZE)]
            for i, line in enumerate(pos_lines): 
                if i % 2 == 0:
                    continue
                donor_motif = line[200:202]
                acceptor_motif = line[598:600]
                add_motif(donor_motif, acceptor_motif, 'pos')

            overall_pos.extend(pos_lines)
    
    with open(pos_file, 'w') as fw: 
        fw.writelines(overall_pos)

# RUNNER
neg_file = f'{OUTPUT_DIR}neg_input.fa'
pos_file = f'{OUTPUT_DIR}pos_input.fa'
sample(NEG_DATASETS, POS_DATASETS, neg_file, pos_file)
display_stats()