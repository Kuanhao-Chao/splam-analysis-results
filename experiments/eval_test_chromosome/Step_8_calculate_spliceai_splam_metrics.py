import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, auc
import seaborn as sns

################################
# SpliceAI output
################################
scores_csv="/ccb/cybertron/khchao/splam-analysis-results/results/eval_test_chromosome/merged.csv"

scores_df = pd.read_csv(scores_csv, sep='\t', header=0)

# pos_scores_df = scores_df[scores_df['ref_label'] == 1]
# print("pos_scores_df: ", len(pos_scores_df))



################################
# Calculating minimum of the scores
################################
scores_df['spliceai_score_min'] = np.minimum(scores_df['spliceai_score_a'], scores_df['spliceai_score_d'])
scores_df['splam_score_min'] = np.minimum(scores_df['splam_score_a'], scores_df['splam_score_d'])

# scores_df['spliceai_score_min'] = (scores_df['spliceai_score_a'] + scores_df['spliceai_score_d']) / 2
# scores_df['splam_score_min'] = (scores_df['splam_score_a'] + scores_df['splam_score_d']) / 2

print("scores_df: ", scores_df)


# Assuming a threshold (you may need to adjust this)
threshold = 0.5

# Define a function to calculate metrics
def calculate_metrics(df, score_col, label_col, threshold):
    # Binarize predictions based on threshold
    predictions = (df[score_col] >= threshold).astype(int)
    
    TP = ((predictions == 1) & (df[label_col] == 1)).sum()
    TN = ((predictions == 0) & (df[label_col] == 0)).sum()
    FP = ((predictions == 1) & (df[label_col] == 0)).sum()
    FN = ((predictions == 0) & (df[label_col] == 1)).sum()
    
    precision = TP / (TP + FP) if TP + FP != 0 else 0
    recall = TP / (TP + FN) if TP + FN != 0 else 0
    accuracy = (TP + TN) / (TP + TN + FP + FN)
    f1_score = 2 * (precision * recall) / (precision + recall) if precision + recall != 0 else 0
    
    return precision, recall, accuracy, f1_score

# Calculate metrics for spliceai_score_min
spliceai_metrics = calculate_metrics(scores_df, 'spliceai_score_min', 'ref_label', threshold)
print(f"SpliceAI metrics (Precision, Recall, Accuracy, F1): {spliceai_metrics}")

# Calculate metrics for splam_score_min
splam_metrics = calculate_metrics(scores_df, 'splam_score_min', 'ref_label', threshold)
print(f"Splam metrics (Precision, Recall, Accuracy, F1): {splam_metrics}")



# ################################
# # Compute ROC and PR curves
# ################################
# fpr_spliceai, tpr_spliceai, _ = roc_curve(scores_df['ref_label'], scores_df['spliceai_score_min'])
# roc_auc_spliceai = auc(fpr_spliceai, tpr_spliceai)

# precision_spliceai, recall_spliceai, _ = precision_recall_curve(scores_df['ref_label'], scores_df['spliceai_score_min'])
# pr_auc_spliceai = auc(recall_spliceai, precision_spliceai)

# fpr_splam, tpr_splam, _ = roc_curve(scores_df['ref_label'], scores_df['splam_score_min'])
# roc_auc_splam = auc(fpr_splam, tpr_splam)

# precision_splam, recall_splam, _ = precision_recall_curve(scores_df['ref_label'], scores_df['splam_score_min'])
# pr_auc_splam = auc(recall_splam, precision_splam)

# ################################
# # Plot ROC curves
# ################################
# plt.figure(figsize=(10, 5))
# plt.subplot(1, 2, 1)
# plt.plot(fpr_spliceai, tpr_spliceai, color='darkorange', lw=2, label='ROC curve (area = %0.2f)' % roc_auc_spliceai)
# plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive Rate')
# plt.title('ROC Curve for spliceai score')
# plt.legend(loc="lower right")

# plt.subplot(1, 2, 2)
# plt.plot(fpr_splam, tpr_splam, color='darkorange', lw=2, label='ROC curve (area = %0.2f)' % roc_auc_splam)
# plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive Rate')
# plt.title('ROC Curve for splam score')
# plt.legend(loc="lower right")

# plt.tight_layout()
# plt.savefig('figure_spliceai_vs_splam_roc.png')

# ################################
# # Plot PR curves
# ################################
# plt.figure(figsize=(10, 5))
# plt.subplot(1, 2, 1)
# plt.plot(recall_spliceai, precision_spliceai, color='darkorange', lw=2, label='PR curve (area = %0.2f)' % pr_auc_spliceai)
# plt.xlabel('Recall')
# plt.ylabel('Precision')
# plt.title('PR Curve for spliceai score')
# plt.legend(loc="lower left")

# plt.subplot(1, 2, 2)
# plt.plot(recall_splam, precision_splam, color='darkorange', lw=2, label='PR curve (area = %0.2f)' % pr_auc_splam)
# plt.xlabel('Recall')
# plt.ylabel('Precision')
# plt.title('PR Curve for splam score')
# plt.legend(loc="lower left")

# plt.tight_layout()
# plt.savefig('figure_spliceai_vs_splam_pr.png')