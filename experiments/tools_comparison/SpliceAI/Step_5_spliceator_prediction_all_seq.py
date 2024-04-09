import tensorflow as tf
# from tensorflow.keras import backend as K
import os, sys
import h5py
import matplotlib.pyplot as plt
from utils import *
import gc

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib.pyplot as plt; plt.rcdefaults()
from tqdm import tqdm
import warnings
import keras
from pkg_resources import resource_filename
from pangolin.model import *

CONSTANT_SIZE = 10
def seq_name(seq):
    return seq[1:]

project_root = '/ccb/cybertron/khchao/splam-analysis-results/' 
project_dir = f'{project_root}/train/results/tool_benchmark/spliceai/'

def main(argv):
    BATCH_SIZE = 100
    BATCH_SIZE_BASE = 1000
    TYPE = "noN"
    all_lines = []
    label = '.'
    da_faf = f'{project_root}train/results/tool_benchmark/spliceai/spliceai.{TYPE}.juncs.seq.fa'

    # Change this to the desired models. The model that each number corresponds to is listed below.
    tissue_id = 6
    model_nums = [tissue_id]
    # 0 = Heart, P(splice)
    # 1 = Heart, usage
    # 2 = Liver, P(splice)
    # 3 = Liver, usage
    # 4 = Brain, P(splice)
    # 5 = Brain, usage
    # 6 = Testis, P(splice)
    # 7 = Testis, usage

    tissue_names = ['Heart', 'Liver', 'Brain', 'Testis']
    tissue_name = tissue_names[tissue_id//2]
    os.makedirs(f'{project_dir}predict_out/pangolin_{tissue_name}/', exist_ok=True)
    # Your existing setup code here...
    combined_output_path = f'{project_dir}predict_out/pangolin_{tissue_name}/spliceai_all_seq.combined.{TYPE}.tsv'
    os.makedirs(os.path.dirname(combined_output_path), exist_ok=True)
    pidx = 0


    # Change this to the desired sequences and strand for each sequence. If the sequence is N bases long, Pangolin will
    # return scores for the middle N-10000 bases (so if you are interested in the score for a single site, the input should
    # be: 5000 bases before the site, base at the site, 5000 bases after the site). Sequences < 10001 bases can be padded with 'N'.
    # seqs = [10001*'A']

    # Load models
    models = []
    for i in model_nums:
        for j in range(1, 6):
            model = Pangolin(L, W, AR)
            if torch.cuda.is_available():
                model.cuda()
                weights = torch.load(resource_filename("pangolin","models/final.%s.%s.3" % (j, i)))
            else:
                weights = torch.load(resource_filename("pangolin","models/final.%s.%s.3" % (j, i)),
                                    map_location=torch.device('cpu'))
            model.load_state_dict(weights)
            model.eval()
            models.append(model)

    # Get scores
    IN_MAP = np.asarray([[0, 0, 0, 0],
                        [1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, 1, 0],
                        [0, 0, 0, 1]])
    INDEX_MAP = {0:1, 1:2, 2:4, 3:5, 4:7, 5:8, 6:10, 7:11}

    def one_hot_encode(seq, strand):
        seq = seq.upper().replace('A', '1').replace('C', '2')
        seq = seq.replace('G', '3').replace('T', '4').replace('N', '0')
        if strand == '+':
            seq = np.asarray(list(map(int, list(seq))))
        elif strand == '-':
            seq = np.asarray(list(map(int, list(seq[::-1]))))
            seq = (5 - seq) % 5  # Reverse complement
        return IN_MAP[seq.astype('int8')]



    

    with open(combined_output_path, "a") as output_fw:
        # Your existing FASTA processing code here...
        #################################
        ## Get all lines of sequences in FASTA file
        #################################
        with open(da_faf, "r") as f:
            print("Processing : ", da_faf)
            lines = f.read().splitlines()
            all_lines = lines
        print("all_lines  : ", len(all_lines))

        #################################
        ## Processing 'POSITIVE' samples
        #################################
        while pidx < len(all_lines):
            if pidx % 2 == 0:
                header_parts = all_lines[pidx].split(";")
                chr = header_parts[0][1:]
                start, end, strand = header_parts[1:]
                junction_name = f'JUNC_{pidx // 2}'
            elif pidx % 2 == 1:
                seq = all_lines[pidx]
                print("seq: ", len(seq))
                # print("Sequence %d: %s" % (i, seq))
                seq = one_hot_encode(seq, strand).T
                seq = torch.from_numpy(np.expand_dims(seq, axis=0)).float()

                if torch.cuda.is_available():
                    seq = seq.to(torch.device("cuda"))

                for j, model_num in enumerate(model_nums):
                    score = []
                    # Average across 5 models
                    for model in models[5*j:5*j+5]:
                        with torch.no_grad():
                            score.append(model(seq)[0][INDEX_MAP[model_num],:].cpu().numpy())
                    mean_score = np.mean(score, axis=0)
                    donor_p = mean_score[199]  # Adjusted indices if necessary
                    acceptor_p = mean_score[-200]
                    print("donor_p: ", donor_p)
                    print("acceptor_p: ", acceptor_p)
                    # Write the combined line in the desired format
                    output_line = f'{chr}\t{start}\t{end}\t{junction_name}\t0.0\t{strand}\t{donor_p}\t{acceptor_p}\n'
                    output_fw.write(output_line)
            pidx += 1

            # if pidx % 4 == 0:
            #     break

if __name__ == "__main__":
    main(sys.argv[1:])
