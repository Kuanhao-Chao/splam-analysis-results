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
project_dir = f'{project_root}results/eval_test_chromosome/spliceai/'



def main(argv):
    BATCH_SIZE = 100
    BATCH_SIZE_BASE = 1000
    # TYPE = argv[1]
    # output_file = argv[2]
    spliceai_version = argv[0]
    print("spliceai_version: ", spliceai_version)
    TYPE = "noN"
    model_path = f'{project_root}experiments/eval_test_chromosome/spliceai_models/spliceai{spliceai_version}.h5'
    print(">> model_path\t\t: ", model_path)
    model = keras.saving.load_model(model_path)
    all_lines = []
    label = '.'
    da_faf = f'{project_root}results/eval_test_chromosome/spliceai/spliceai.{TYPE}.juncs.seq.fa'

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



# def main(argv):
#     BATCH_SIZE = 100
#     BATCH_SIZE_BASE = 1000
#     # TYPE = argv[1]
#     # output_file = argv[2]
#     # spliceai_version = argv[3]
#     TYPE = "noN"
#     spliceai_version = 2
#     model_path = f'{project_root}experiments/eval_test_chromosome/spliceai_models/spliceai{spliceai_version}.h5'
#     print(">> model_path\t\t: ", model_path)
#     model = keras.saving.load_model(model_path)
#     all_lines = []
#     label = '.'
#     da_faf = f'{project_root}results/eval_test_chromosome/spliceai/spliceai.{TYPE}.juncs.seq.fa'

#     os.makedirs(f'{project_dir}predict_out/spliceai{spliceai_version}/', exist_ok=True)
#     d_score_tsv_f = f'{project_dir}predict_out/spliceai{spliceai_version}/spliceai_all_seq.score.d.{TYPE}.tsv'
#     a_score_tsv_f = f'{project_dir}predict_out/spliceai{spliceai_version}/spliceai_all_seq.score.a.{TYPE}.tsv'
#     n_score_tsv_f = f'{project_dir}predict_out/spliceai{spliceai_version}/spliceai_all_seq.score.n.{TYPE}.tsv'
#     name_tsv_f = f'{project_dir}predict_out/spliceai{spliceai_version}/spliceai_all_seq.name.{TYPE}.tsv'

#     d_score_fw = open(d_score_tsv_f, "a")
#     a_score_fw = open(a_score_tsv_f, "a")
#     n_score_fw = open(n_score_tsv_f, "a")
#     name_fw = open(name_tsv_f, "a")

#     #################################
#     ## Get all lines of sequences in FASTA file
#     #################################
#     with open(da_faf, "r") as f:
#         print("Processing : ", da_faf)
#         lines = f.read().splitlines()
#         all_lines = lines
#     print("all_lines  : ", len(all_lines))

#     #################################
#     ## Processing 'POSITIVE' samples
#     #################################
#     COUNTER = 0
#     pidx = 0
#     pidx = BATCH_SIZE-BATCH_SIZE_BASE
#     seq = ""
#     print("BATCH_SIZE     : ", BATCH_SIZE)
#     print("BATCH_SIZE_BASE: ", BATCH_SIZE_BASE)
#     while pidx < len(all_lines):
#         if pidx % 2 == 0:
#             chr, start, end, strand = all_lines[pidx].split(";")
#             chr = chr[1:]
#             name_fw.write(' '.join([chr, start, end, strand])+"\n")
#             pass
#         elif pidx % 2 == 1:
#             seq = all_lines[pidx]
#             X, Y = create_datapoints(seq, label)
#             X = X[None, :]
#             Y = np.array(Y)
#             X = tf.convert_to_tensor(X, dtype=tf.float32)
#             Y = tf.convert_to_tensor(Y, dtype=tf.float32)
#             Y_pred = model.predict(X)
#             COUNTER += 1
#             print("X.shape     : ", X.shape)
#             print("Y_pred.shape: ", Y_pred.shape)
#             donor_p = Y_pred[0][200-1][2]
#             acceptor_p = Y_pred[0][len(Y_pred)-200-1][1]
#             print("(chr, start, end, strand): ", (chr, start, end, strand))
#             print("donor_p    : ", donor_p)
#             print("acceptor_p : ", acceptor_p)
#             Y_pred = np.array(Y_pred)
#             print(Y_pred)
#             print("Y_pred.shape: ", Y_pred.shape)
#             d_score_fw.write(str(donor_p)+"\n")
#             a_score_fw.write(str(acceptor_p)+"\n")            
#         pidx += 1
#         print("====================")
#         if pidx %100 == 0:
#             print("pidx: ", pidx)
#     d_score_fw.close()
#     a_score_fw.close()
#     n_score_fw.close()
#     name_fw.close()

if __name__ == "__main__":
    main(sys.argv[1:])
