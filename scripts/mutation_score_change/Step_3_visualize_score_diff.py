import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
import sys, os
import pandas as pd
import numpy as np

type = sys.argv[1]
target = "pos_MANE"
project_root = "../../train/src/INPUTS/800bp/"
output_dir = f"../../results/mutation/{target}/scores/"

def plot_violin(df, target):
    plt.figure(figsize=(80, 6))
    sns.violinplot(x='position', y='diff', data=df)
    plt.xlabel('Position')
    plt.ylabel(f'{target} score difference')
    plt.title(f'Violin Plot of score difference for {target} Sites at each Position')
    plt.savefig(f'{output_dir}{target}score_difference_violin.png')


def plot_average(df, target):
    plt.figure(figsize=(15, 6))
    sns.lineplot(x='position', y='diff', data=df)
    plt.xlabel('Position')
    plt.ylabel(f'Average {target} score difference')
    plt.title(f'Average Score Difference for {target} Sites at Each Position')
    plt.savefig(f'{output_dir}{target}_score_difference_average.png')


def main(fasta_file):
    # Create a DataFrame to store d_diff for each position
    all_d_diffs = pd.DataFrame()
    all_a_diffs = pd.DataFrame()

    org_d_a_score_bed = f"{output_dir}original/LOG/_junction_score.bed"
    org_df = pd.read_csv(org_d_a_score_bed, sep="\t", header=None, names=["chr", "start", "end", "name", "score", "strand", "donor", "acceptor"])
    # print(org_df)

    for position in range(800):
        d_a_score_bed = f"{output_dir}pos_{position}/LOG/_junction_score.bed"
        df = pd.read_csv(d_a_score_bed, sep="\t", header=None, names=["chr", "start", "end", "name", "score", "strand", "donor", "acceptor"])
        # print(df)

        # Joining the DataFrames
        result = pd.merge(org_df, df, on=['chr', 'start', 'end', 'strand'])
        print(result)

        d_diff = result["donor_y"] - result["donor_x"]
        a_diff = result["acceptor_y"] - result["acceptor_x"]
        # Add position and d_diff to all_d_diffs DataFrame
        temp_d_df = pd.DataFrame({'position': position, 'diff': d_diff})
        all_d_diffs = pd.concat([all_d_diffs, temp_d_df])
        
        # Add position and a_diff to all_a_diffs DataFrame
        temp_a_df = pd.DataFrame({'position': position, 'diff': a_diff})
        all_a_diffs = pd.concat([all_a_diffs, temp_a_df])
        # if position == 110:
        #     break
    # print(all_d_diffs)
    # print(all_a_diffs)
    # plot_violin(all_d_diffs, "Donor")
    # plot_violin(all_a_diffs, "Acceptor")

    # Calculate average differences
    avg_d_diffs = all_d_diffs.groupby('position').mean().reset_index()
    avg_a_diffs = all_a_diffs.groupby('position').mean().reset_index()

    # print(avg_d_diffs["diff"])
    # print(avg_a_diffs["diff"][avg_a_diffs["diff"] < -0.01])

    plot_average(avg_d_diffs, "Donor")
    plot_average(avg_a_diffs, "Acceptor")



if __name__ == "__main__":
    os.makedirs(output_dir, exist_ok=True)
    fasta_file = f"{project_root}input_{target}/{type}_{target}.shuffle.fa"  # Replace with your FASTA file path
    main(fasta_file)