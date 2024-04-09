from Bio import SeqIO
import random
import os
import sys
from pathlib import Path
import pandas as pd

project_root = '/ccb/cybertron/khchao/splam-analysis-results'
def main():
    SEQ_LEN = 800
    THRESHOLD = 100
    ref_junc_f = f'{project_root}/train/results/ALL_RefSeq/REF_junctions/ref_d_a.sort.bed'
    test_junc_f = f'{project_root}/train/results/tool_benchmark/{SEQ_LEN}bp/d_a.bed'
    junc_fs = [test_junc_f]
    print('junc_fs: ', junc_fs)

    output_dir = f'{project_root}/train/results/tool_benchmark/'
    output_files = [output_dir]
    print('output_files: ', output_files)
    COUNTER = 0
    column_name = ['chr', 'start', 'end', 'junc', 'score', 'strand', 'label']
    for junc_fidx in range(0, len(junc_fs), 1):
        print("junc_fidx\t: ", junc_fidx)
        junc_f = junc_fs[junc_fidx]
        print("junc_f\t\t: ", junc_f)
        os.makedirs(output_files[junc_fidx]+"splam/", exist_ok=True)
        os.makedirs(output_files[junc_fidx]+"spliceai/", exist_ok=True)
        junc_df = pd.read_csv(junc_f, delimiter="\t", header=None)
        junc_df = junc_df.loc[junc_df[1] > 0]

        ################################
        # SPLAM test data curatiolsn
        ################################
        junc_df.columns = column_name
        # junc_df_spliceai = junc_df.copy()
        print('junc_df: ', junc_df)
        print('len(junc_df): ', len(junc_df))
        junc_df['end'] -= 1
        print("SPLAM junc_df   : ", junc_df)
        junc_df.to_csv(output_files[junc_fidx]+"splam/splam.juncs.bed", sep="\t", header=None, index=0)

        ################################
        # SpliceAI test data curatiolsn
        ################################
        junc_df['start'] -= 5200
        junc_df['end'] += 5200
        print("SpliceAI out_df: ", junc_df)
        junc_df.to_csv(output_files[junc_fidx]+"spliceai/spliceai.juncs.bed", sep="\t", header=None, index=0)

if __name__ == "__main__":
    main()