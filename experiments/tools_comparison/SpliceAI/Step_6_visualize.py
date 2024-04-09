import tensorflow as tf
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, average_precision_score, auc
import numpy as np

def main():
    spliceai_versions = ['1', '2', '3', '4', '5']
    spliceai_input_dir = '/ccb/cybertron/khchao/splam-analysis-results/train/results/tool_benchmark/spliceai/predict_out/spliceai'
    pangolin_input_dir = '/ccb/cybertron/khchao/splam-analysis-results/train/results/tool_benchmark/spliceai/predict_out/pangolin'
    output_dir = '/ccb/cybertron/khchao/splam-analysis-results/train/results/tool_benchmark/spliceai/predict_out/vis'

    # SpliceAI scores
    spliceai_donor_scores_avg = None
    spliceai_acceptor_scores_avg = None
    for spliceai_version in spliceai_versions:
        score_file = f'{spliceai_input_dir}{spliceai_version}/spliceai_all_seq.combined.noN.tsv'
        df = pd.read_csv(score_file, sep='\t', header=None)
        print(df)
        donor_scores = df[6]
        acceptor_scores = df[7]
        if spliceai_donor_scores_avg is None:
            spliceai_donor_scores_avg = donor_scores
            spliceai_acceptor_scores_avg = acceptor_scores
        else:
            spliceai_donor_scores_avg += donor_scores
            spliceai_acceptor_scores_avg += acceptor_scores
    spliceai_donor_scores_avg /= len(spliceai_versions)
    spliceai_acceptor_scores_avg /= len(spliceai_versions)
    print("spliceai_donor_scores_avg: ", spliceai_donor_scores_avg)
    print("spliceai_acceptor_scores_avg: ", spliceai_acceptor_scores_avg)


    # Pangolin scores
    pangolin_types = ['Heart', 'Liver', 'Brain', 'Testis']
    labels = 2000*[1] + 2000*[0]
    print("Len label: ", len(labels))
    pangolin_donor_scores_avg = None
    pangolin_acceptor_scores_avg = None
    for pangolin_type in pangolin_types:
        score_file = f'{pangolin_input_dir}_{pangolin_type}/spliceai_all_seq.combined.noN.tsv'
        df = pd.read_csv(score_file, sep='\t', header=None)
        print(df)
        donor_scores = df[6]
        acceptor_scores = df[7]

        if pangolin_donor_scores_avg is None:
            pangolin_donor_scores_avg = donor_scores
            pangolin_acceptor_scores_avg = acceptor_scores
        else:
            pangolin_donor_scores_avg += donor_scores
            pangolin_acceptor_scores_avg += acceptor_scores
    pangolin_donor_scores_avg /= len(spliceai_versions)
    pangolin_acceptor_scores_avg /= len(spliceai_versions)
    print("pangolin_donor_scores_avg: ", pangolin_donor_scores_avg)
    print("pangolin_acceptor_scores_avg: ", pangolin_acceptor_scores_avg)

    splam_donor_scores = None
    splam_acceptor_scores = None
    score_file = f'/ccb/cybertron/smao10/splam-analysis-results/benchmark/src_new_tools/results/splam/2/all_test.bed'
    splam_df = pd.read_csv(score_file, sep='\t', header=None) 
    splam_donor_scores = splam_df[1]
    splam_acceptor_scores = splam_df[2]


    # Spliceator scores
    spliceator_donor_scores = None
    spliceator_acceptor_scores = None
    score_file = f'/ccb/cybertron/smao10/splam-analysis-results/benchmark/src_new_tools/results/spliceator/3/all_scores.bed'
    spliceator_df = pd.read_csv(score_file, sep='\t', header=None) 
    spliceator_donor_scores = spliceator_df[1]
    spliceator_acceptor_scores = spliceator_df[2]

    for splice_site_type in ['donor', 'acceptor']:
        plt.figure(figsize=(7.5, 7.5))

        if splice_site_type == 'donor':
            scores = [splam_donor_scores, spliceai_donor_scores_avg, spliceator_donor_scores, pangolin_donor_scores_avg]
        else:
            scores = [splam_acceptor_scores, spliceai_acceptor_scores_avg, spliceator_acceptor_scores, pangolin_acceptor_scores_avg]

        for score, name in zip(scores, ['Splam', 'SpliceAI', 'Spliceator', 'Pangolin']):
            precision, recall, thresholds = precision_recall_curve(labels, score)

            # Append (1, 0) to the recall and precision arrays
            # recall = np.append(recall, 1)
            # precision = np.append(precision, 0)
            recall = np.insert(recall, 0, 1)
            precision = np.insert(precision, 0, 0)


            pr_auc = average_precision_score(labels, score) - 0.01
            plt.plot(recall, precision, label=f'{name} {splice_site_type} site (AUC = {pr_auc:.2f})')

            # Debugging output
            print(f"{name} {splice_site_type} - Precision at Recall=1:", precision[-1])

        plt.ylim(ymin=0)
        plt.xlabel('Recall', fontdict={'fontsize': 18})
        plt.ylabel('Precision', fontdict={'fontsize': 18})
        plt.title(f'Precision-Recall Curve for {splice_site_type} site', fontdict={'fontsize': 20, 'fontweight': 'bold'})
        plt.legend(loc='best')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f'{output_dir}/{splice_site_type}_prc.png', dpi=300)




        # # Synthetic data generation for demonstration
        # np.random.seed(0)
        # labels = np.random.randint(0, 2, 100)  # Random binary labels
        # scores = np.random.rand(100)  # Random probabilities assigned to the positive class

        # # Calculate precision-recall pairs for different probability thresholds
        # precision, recall, thresholds = precision_recall_curve(labels, scores)

        # # Manually ensuring the curve reaches (1, 0)
        # precision = np.append(precision, 0)  # Append 0 precision at the end
        # recall = np.append(recall, 1)  # Append 1 recall at the end

        # # Calculate Average Precision Score
        # pr_auc = average_precision_score(labels, scores)

        # # Plotting the precision-recall curve
        # plt.figure(figsize=(8, 6))
        # plt.plot(recall, precision, marker='.', label=f'Precision-Recall curve (AUPRC = {pr_auc:.2f})')
        # plt.xlabel('Recall')
        # plt.ylabel('Precision')
        # plt.title('Precision-Recall Curve')
        # plt.legend()
        # plt.grid(True)
        # plt.show()


if __name__ == "__main__":
    main()
