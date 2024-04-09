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

# Import
from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice import MMSplice, predict_save, predict_all_table
from mmsplice.utils import max_varEff
from mmsplice.layers import ConvDNA

# print("ConvDNA: ", ConvDNA)

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



    # example files
    repo_root_dir = '/ccb/cybertron/khchao/splam-analysis-results/experiments/tools_comparison/MMSplice/MMSplice_MTSplice/'
    gtf = f'{repo_root_dir}tests/data/test.gtf'
    vcf = f'{repo_root_dir}tests/data/test.vcf.gz'
    fasta = f'{repo_root_dir}tests/data/hg19.nochr.chr17.fa'
    csv = 'pred.csv'
    print(">> gtf\t\t: ", gtf)

    dl = SplicingVCFDataloader(gtf, fasta, vcf, tissue_specific=False)
    # Specify model
    print("Hello world!")
    print("MMSplice: ", MMSplice)
    model = MMSplice()
    # model = MMSplice(
    #     **{k: v for k, v in options.items() if v})

    # # Or predict and return as df
    # predictions = predict_all_table(model, dl, pathogenicity=True, splicing_efficiency=True)

if __name__ == "__main__":
    main(sys.argv[1:])
