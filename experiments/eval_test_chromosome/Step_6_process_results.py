import pandas as pd

project_root="/ccb/cybertron/khchao/splam-analysis-results/"

################################
# SpliceAI output
################################
spliceai_f = "/ccb/cybertron/khchao/splam-analysis-results/results/eval_test_chromosome/spliceai/predict_out/spliceai3/spliceai_all_seq.combined.noN.tsv"
spliceai_df = pd.read_csv(spliceai_f, sep='\t', header=None)
spliceai_df.columns = ['chromosome', 'start', 'end', 'name', 'label', 'strand', 'score_a', 'score_d']
print("spliceai_df: ", spliceai_df)

################################
# Splam output
################################
splam_f = "/ccb/cybertron/khchao/splam-analysis-results/results/eval_test_chromosome/splam/predict_out/splam/LOG/_junction_score.bed"
splam_df = pd.read_csv(splam_f, sep='\t', header=None)
splam_df.columns = ['chromosome', 'start', 'end', 'name', 'label', 'strand', 'score_a', 'score_d']
print("splam_df: ", splam_df)

################################
# Merged dataframe
################################
# Merge the two dataframes by ['chromosome', 'start', 'end', 'strand']
spliceai_splam_df = pd.merge(spliceai_df, splam_df, on=['chromosome', 'start', 'end', 'strand'], how='inner')
spliceai_splam_df.columns = ['chromosome', 'start', 'end', 'name_spliceai', 'label_x', 'strand', 'spliceai_score_a', 'spliceai_score_d', 'name_splam', 'label_y', 'splam_score_a', 'splam_score_d']
# Display the merged DataFrame
merged_df = spliceai_splam_df.drop(['name_splam', 'label_x', 'name_spliceai', 'label_y'], axis=1)
print(merged_df)

################################
# Get label
################################
ref_junc_f = f'{project_root}/train/results/ALL_RefSeq/REF_junctions/ref_d_a_all.sort.bed'
ref_junc_df = pd.read_csv(ref_junc_f, sep='\t', header=None)
ref_junc_df.columns = ['chromosome', 'start', 'end', 'name', 'ref_label', 'strand']
# ref_junc_df['start'] = ref_junc_df['start'] + 1
ref_junc_df['end'] = ref_junc_df['end'] - 1
print("ref_junc_df: ", ref_junc_df)
final_df = pd.merge(merged_df, ref_junc_df, on=['chromosome', 'start', 'end', 'strand'], how='left')
final_df.drop(columns=['name'], inplace=True)
final_df['ref_label'] = final_df['ref_label'].fillna(0)
# Display the merged DataFrame
print("final_df: ", final_df)

print(len(final_df[final_df['ref_label'] == 1]))

final_df.to_csv("/ccb/cybertron/khchao/splam-analysis-results/results/eval_test_chromosome/merged.csv", sep='\t', index=False)


################################
# Score visualization
################################
import matplotlib.pyplot as plt
import seaborn as sns
# Create scatter plots for donor and acceptor scores comparison
plt.figure(figsize=(12, 6))

# Scatter plot for donor scores comparison
plt.subplot(1, 2, 1)
plt.scatter(final_df['spliceai_score_a'], final_df['splam_score_a'], c='red', alpha=0.4)
plt.scatter(final_df.loc[final_df['ref_label'] == 1, 'spliceai_score_a'], 
            final_df.loc[final_df['ref_label'] == 1, 'splam_score_a'], 
            c='darkgreen', alpha=0.8)
plt.title('Donor Scores Comparison')
plt.xlabel('SpliceAI Donor Score')
plt.ylabel('Splam Donor Score')

# Scatter plot for acceptor scores comparison
plt.subplot(1, 2, 2)
plt.scatter(final_df['spliceai_score_d'], final_df['splam_score_d'], c='red', alpha=0.4)
plt.scatter(final_df.loc[final_df['ref_label'] == 1, 'spliceai_score_d'], 
            final_df.loc[final_df['ref_label'] == 1, 'splam_score_d'], 
            c='darkgreen', alpha=0.8)
plt.title('Acceptor Scores Comparison')
plt.xlabel('SpliceAI Acceptor Score')
plt.ylabel('Splam Acceptor Score')

# Show plots
plt.tight_layout()
plt.savefig('/ccb/cybertron/khchao/splam-analysis-results/results/eval_test_chromosome/figure_spliceai_vs_splam_scores_colored.png', dpi=300)
plt.clf()

# Create scatter plots for SpliceAI and Splam scores comparison
plt.figure(figsize=(12, 6))

# Scatter plot for SpliceAI scores comparison
plt.subplot(1, 2, 1)
plt.scatter(final_df['spliceai_score_d'], final_df['spliceai_score_a'], c='red', alpha=0.4)
plt.scatter(final_df.loc[final_df['ref_label'] == 1, 'spliceai_score_d'], 
            final_df.loc[final_df['ref_label'] == 1, 'spliceai_score_a'], 
            c='darkgreen', alpha=0.8)
plt.title('SpliceAI Scores Comparison')
plt.xlabel('SpliceAI Donor Score')
plt.ylabel('SpliceAI Acceptor Score')

# Scatter plot for Splam scores comparison
plt.subplot(1, 2, 2)
plt.scatter(final_df['splam_score_d'], final_df['splam_score_a'], c='red', alpha=0.4)
plt.scatter(final_df.loc[final_df['ref_label'] == 1, 'splam_score_d'], 
            final_df.loc[final_df['ref_label'] == 1, 'splam_score_a'], 
            c='darkgreen', alpha=0.8)
plt.title('Splam Scores Comparison')
plt.xlabel('Splam Donor Score')
plt.ylabel('Splam Acceptor Score')

# Show plots
plt.tight_layout()
plt.savefig('/ccb/cybertron/khchao/splam-analysis-results/results/eval_test_chromosome/figure_donor_vs_acceptor_colored.png', dpi=300)
plt.clf()