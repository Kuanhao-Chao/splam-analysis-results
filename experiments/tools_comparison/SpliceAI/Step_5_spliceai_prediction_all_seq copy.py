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

CONSTANT_SIZE = 10
def seq_name(seq):
    return seq[1:]

project_root = '/ccb/cybertron/khchao/splam-analysis-results/' 
project_dir = f'{project_root}/train/results/tool_benchmark/spliceai/'

def main(argv):
    BATCH_SIZE = 100
    BATCH_SIZE_BASE = 1000
    spliceai_version = argv[0]
    print("spliceai_version: ", spliceai_version)
    TYPE = "noN"
    model_path = f'{project_root}experiments/eval_test_chromosome/spliceai_models/spliceai{spliceai_version}.h5'
    print(">> model_path\t\t: ", model_path)
    model = keras.saving.load_model(model_path)
    all_lines = []
    label = '.'
    da_faf = f'{project_root}train/results/tool_benchmark/spliceai/spliceai.{TYPE}.juncs.seq.fa'
    os.makedirs(f'{project_dir}predict_out/spliceai{spliceai_version}/', exist_ok=True)
    # Your existing setup code here...
    model_path = f'{project_root}experiments/eval_test_chromosome/spliceai_models/spliceai{spliceai_version}.h5'    
    combined_output_path = f'{project_dir}predict_out/spliceai{spliceai_version}/spliceai_all_seq.combined.{TYPE}.tsv'
    os.makedirs(os.path.dirname(combined_output_path), exist_ok=True)
    pidx = 0

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
                X, _ = create_datapoints(seq, label)
                X = tf.expand_dims(X, axis=0)
                Y_pred = model.predict(X)
                donor_p = Y_pred[0][199][2]  # Adjusted indices if necessary
                acceptor_p = Y_pred[0][-200][1]
                # Write the combined line in the desired format
                output_line = f'{chr}\t{start}\t{end}\t{junction_name}\t0.0\t{strand}\t{donor_p}\t{acceptor_p}\n'
                output_fw.write(output_line)
            pidx += 1



if __name__ == "__main__":
    main(sys.argv[1:])
