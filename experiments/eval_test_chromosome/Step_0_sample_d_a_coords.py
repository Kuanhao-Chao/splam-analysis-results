from Bio import SeqIO
import random
import os
import sys
from pathlib import Path
import pandas as pd

project_root = '/ccb/cybertron/khchao/splam-analysis-results/'
def main():
    SEQ_LEN = 800
    THRESHOLD = 100
    test_junc_f = f'{project_root}/train/results/BAM_TEST_junctions/{SEQ_LEN}bp/d_a.bed'
    test_junc_out_f = f'{project_root}/train/results/BAM_TEST_junctions/{SEQ_LEN}bp/d_a.sample1.bed'

    # Read the BED file
    df = pd.read_csv(test_junc_f, sep='\t', header=None)
    # df_sampled = df.sample(n=20000, replace=False) if len(df) >= 20000 else df
    df_sampled = df.sample(n=200000, replace=False)
    # Write the sampled entries to a new BED file
    df_sampled.to_csv(test_junc_out_f, sep='\t', index=False, header=False)

if __name__ == "__main__":
    main()