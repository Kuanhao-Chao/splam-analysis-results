import matplotlib.pyplot as plt
import pickle
import numpy as np
import os
import pandas as pd
from util import *
from sklearn.metrics import auc, accuracy_score, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve, PrecisionRecallDisplay
from sklearn import svm


# define a function to calculate precision, recall, F1 score, and queue rate
def calculate_metrics(threshold, true_labels, predict_probabilities):
    predictions = (predict_probabilities >= threshold).astype(int)
    true_positives = np.sum((predictions == 1) & (true_labels == 1))
    false_positives = np.sum((predictions == 1) & (true_labels == 0))
    false_negatives = np.sum((predictions == 0) & (true_labels == 1))
    true_negatives = np.sum((predictions == 0) & (true_labels == 0))
    precision = true_positives / (true_positives + false_positives)
    recall = true_positives / (true_positives + false_negatives)
    f1_score = 2 * (precision * recall) / (precision + recall)
    queue_rate = (true_positives + false_positives) / (true_positives + false_positives + true_negatives + false_negatives)
    return precision, recall, f1_score, queue_rate



def main(db):

    #####################################
    # Declaring parameters for probability & prediction array
    #####################################

    POS_NUM = 2500
    NEG_NUM = 25000

    ### noN ###
    # positive
    noN_pos_df = pd.read_csv(f'../../pos_test/6_output/combine/averaged_data.noN.{db}.csv')
    noN_pos_df['label'] = 1
    # negative
    noN_neg_df = pd.read_csv(f'../../neg_test/6_output/combine/averaged_data.noN.{db}.csv')
    noN_neg_df['label'] = 0
    
    noN_pos_df = noN_pos_df.sample(n=POS_NUM, random_state=1)
    noN_neg_df = noN_neg_df.sample(n=NEG_NUM, random_state=1)

    noN_merge_df = pd.concat([noN_pos_df, noN_neg_df], axis=0)


    print("noN_pos_df: ", len(noN_pos_df))
    print("noN_neg_df: ", len(noN_neg_df))
    print("noN_merge_df: ", len(noN_merge_df))
    
    ### N ###
    # positive
    N_pos_df = pd.read_csv(f'../../pos_test/6_output/combine/averaged_data.N.{db}.csv')
    N_pos_df['label'] = 1
    # negative
    N_neg_df = pd.read_csv(f'../../neg_test/6_output/combine/averaged_data.N.{db}.csv')
    N_neg_df['label'] = 0
    
    N_pos_df = N_pos_df.sample(n=POS_NUM, random_state=1)
    N_neg_df = N_neg_df.sample(n=NEG_NUM, random_state=1)

    N_merge_df = pd.concat([N_pos_df, N_neg_df], axis=0)
    
    # load labels and predictions from dataframe
    spliceai_noN_d_label = noN_merge_df['label']
    spliceai_noN_a_label = noN_merge_df['label']
    spliceai_noN_d_pred = noN_merge_df['d_score_spliceai']
    spliceai_noN_a_pred = noN_merge_df['a_score_spliceai']

    spliceai_N_d_label = N_merge_df['label']
    spliceai_N_a_label = N_merge_df['label']
    spliceai_N_d_pred = N_merge_df['d_score_spliceai']
    spliceai_N_a_pred = N_merge_df['a_score_spliceai']

    splam_d_label = N_merge_df['label']
    splam_a_label = N_merge_df['label']
    splam_d_pred = N_merge_df['d_score_splam']
    splam_a_pred = N_merge_df['a_score_splam']

    print(f'SpliceAI noN_d labels:\n\t1:{len(spliceai_noN_d_label[spliceai_noN_d_label==1])}\n\t0:{len(spliceai_noN_d_label[spliceai_noN_d_label==0])}')
    print(f'SpliceAI noN_a labels:\n\t1:{len(spliceai_noN_a_label[spliceai_noN_a_label==1])}\n\t0:{len(spliceai_noN_a_label[spliceai_noN_a_label==0])}')
    print(f'SpliceAI N_d labels:\n\t1:{len(spliceai_N_d_label[spliceai_N_d_label==1])}\n\t0:{len(spliceai_N_d_label[spliceai_N_d_label==0])}')
    print(f'SpliceAI N_a labels:\n\t1:{len(spliceai_N_a_label[spliceai_N_a_label==1])}\n\t0:{len(spliceai_N_a_label[spliceai_N_a_label==0])}')
    print(f'SPLAM d labels:\n\t1:{len(splam_d_label[splam_d_label==1])}\n\t0:{len(splam_d_label[splam_d_label==0])}')
    print(f'SPLAM a labels:\n\t1:{len(splam_a_label[splam_a_label==1])}\n\t0:{len(splam_a_label[splam_a_label==0])}')

    print(f'SpliceAI noN_d pred {len(spliceai_noN_d_pred)}')
    print(f'SpliceAI noN_a pred {len(spliceai_noN_a_pred)}')
    print(f'SpliceAI N_d pred {len(spliceai_N_d_pred)}')
    print(f'SpliceAI N_d pred {len(spliceai_N_a_pred)}')
    print(f'SPLAM d pred {len(splam_d_pred)}')
    print(f'SPLAM a pred {len(splam_a_pred)}')  



    #####################################
    # DT plot visualization
    #####################################

    # generate some sample data for predict probabilities and true labels
    path = f'./{db}_DT_plot_spliceai_noN_{POS_NUM}-{NEG_NUM}.png'
    predict_probabilities = np.minimum(spliceai_noN_d_pred, spliceai_noN_a_pred)
    true_labels = spliceai_noN_a_label
    plot_DT_plot(true_labels, predict_probabilities, path)

    # generate some sample data for predict probabilities and true labels
    path = f'./{db}_DT_plot_splam_{POS_NUM}-{NEG_NUM}.png'
    predict_probabilities = np.minimum(splam_d_pred,splam_a_pred)
    true_labels = splam_d_label
    plot_DT_plot(true_labels, predict_probabilities, path)



def plot_DT_plot(true_labels, predict_probabilities, path):
    # define the range of threshold values to plot
    thresholds = np.arange(0, 0.99, 0.0001)

    # calculate the metrics for each threshold value
    precisions = []
    recalls = []
    f1_scores = []
    # queue_rates = []

    for threshold in thresholds:
        precision, recall, f1_score, queue_rate = calculate_metrics(threshold, true_labels, predict_probabilities)
        precisions.append(precision)
        recalls.append(recall)
        f1_scores.append(f1_score)
        # queue_rates.append(queue_rate)

    print(f1_scores)
    # find the optimal threshold based on F1 score
    optimal_index = np.argmax(f1_scores)
    optimal_threshold = thresholds[optimal_index]
    optimal_precision = precisions[optimal_index]
    optimal_recall = recalls[optimal_index]
    optimal_f1_score = f1_scores[optimal_index]
    # optimal_queue_rate = queue_rates[optimal_index]
    print(f'Optimal threshold: {optimal_threshold}')
    print(f'Optimal_threshold: {optimal_threshold:.4f}, Precision: {optimal_precision:.2f}, Recall: {optimal_recall:.2f}, F1 score: {optimal_f1_score:.2f}')

    # plot the DT plot
    plt.figure(figsize=(6, 3.5))
    plt.plot(thresholds, precisions, label='Precision', linewidth=2)
    plt.plot(thresholds, recalls, label='Recall', linewidth=2)
    plt.plot(thresholds, f1_scores, label='F1 Score', linewidth=2)
    # plt.plot(thresholds, queue_rates, label='Queue Rate')
    plt.axvline(x=optimal_threshold, linestyle='--', color='r', label='Optimal Threshold (maximum F1 score)')
    plt.xlabel('Threshold', size = 13)
    plt.ylabel('Performance Metrics', size = 13)
    plt.xlim(-0.02, 1.02)
    plt.ylim(0.0, 1.04)
    
    plt.legend(loc='upper right', ncol=4, bbox_to_anchor=(0.95, 1.22), prop={'size':7})

    plt.title(f'Optimal_threshold: {optimal_threshold:.4f}, Precision: {optimal_precision:.2f}, Recall: {optimal_recall:.2f}, F1 score: {optimal_f1_score:.2f}', fontsize=9)
    plt.xticks(size = 10)
    plt.yticks(size = 10)
    plt.tight_layout()

    os.makedirs(os.path.dirname(path), exist_ok=True)
    plt.savefig(path, dpi=300)
    plt.close()

if __name__ == "__main__":
    
    databases = ['GRCm39', 'Mmul_10', 'NHGRI_mPanTro3', 'TAIR10']
    idxs = [0,2,3]
    
    for idx in idxs:
        main(databases[idx])
