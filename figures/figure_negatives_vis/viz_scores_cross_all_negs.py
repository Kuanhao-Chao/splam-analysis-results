import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
# Load the data into a DataFrame
from io import StringIO
# Assuming the data is loaded into a DataFrame named `df`

# Initialize lists to hold aggregated data
thresholds_list = list(range(2, 101))
mean_scores_a = []
std_devs_a = []
mean_scores_d = []
std_devs_d = []


# Loop through each threshold to load data, and compute mean and standard deviation
fig, axes = plt.subplots(1, 2, figsize=(9, 6))  # Adjusting the size for better visibility
for threshold in thresholds_list:
    neg_pred_file = f'/ccb/cybertron/khchao/splam-analysis-results/src_tools_evaluation/splam_result/RELEASED/neg_{threshold}/LOG/_junction_score.bed'
    df = pd.read_csv(neg_pred_file, sep="\t", header=None, names=["chr", "start", "end", "junction", "unknown", "strand", "donor_score", "acceptor_score"])
    df = df.sample(n=600)
    # Assuming you're interested in donor_score; replace with acceptor_score if necessary
    filtered_scores_a = df['acceptor_score'][:600]
    mean_scores_a.append(filtered_scores_a.mean())
    std_devs_a.append(filtered_scores_a.std())

    filtered_scores_d = df['donor_score'][:600]
    mean_scores_d.append(filtered_scores_d.mean())
    std_devs_d.append(filtered_scores_d.std())

axes[0].set_title('Donor Site')
axes[0].set_xlabel('Spliced Alignment Count')
axes[0].set_ylabel('Average Splam Score')
# axes[0].errorbar(thresholds_list, mean_scores_d, yerr=std_devs_d, fmt='o', ecolor='lightgray', elinewidth=2, capsize=0, label='Mean Donor Score with Error Bar', markersize = 3)
axes[0].plot(thresholds_list, mean_scores_d, 'o', label='Mean Donor Score', markersize = 3)
# Calculating the trend line
z = np.polyfit(thresholds_list, mean_scores_d, 1)  # Fit a 1st degree polynomial to the data
p = np.poly1d(z)  # Create the polynomial equation based on the fit
# Plotting the trend line
axes[0].plot(thresholds_list, p(thresholds_list), "r--", label='Donor trend line', linewidth=3.0)  # "r--" makes the line red and dashed
axes[0].legend()

axes[1].set_title('Acceptor Site')
axes[1].set_xlabel('Spliced Alignment Counts')
axes[1].set_ylabel('Average Splam Score')
# axes[1].errorbar(thresholds_list, mean_scores_a, yerr=std_devs_a, fmt='o', ecolor='lightgray', elinewidth=2, capsize=0, label='Mean Acceptor Score with Error Bar', markersize = 3)
axes[1].plot(thresholds_list, mean_scores_a, 'o', label='Mean Acceptor Score', markersize = 3)
# Calculating the trend line
z = np.polyfit(thresholds_list, mean_scores_a, 1)  # Fit a 1st degree polynomial to the data
p = np.poly1d(z)  # Create the polynomial equation based on the fit
# Plotting the trend line
axes[1].plot(thresholds_list, p(thresholds_list), "r--", label='Acceptor trend line', linewidth=3.0)  # "r--" makes the 
axes[1].legend()

plt.legend()
plt.tight_layout()
plt.savefig('thresholds.png', dpi=300)
plt.clf()



# Loop through each threshold to load data, and compute mean and standard deviation
filtered_scores_a = []
filtered_scores_d = []
selected_threshold = 0.5
for threshold in thresholds_list:
    neg_pred_file = f'/ccb/cybertron/khchao/splam-analysis-results/src_tools_evaluation/splam_result/RELEASED/neg_{threshold}/LOG/_junction_score.bed'
    df = pd.read_csv(neg_pred_file, sep="\t", header=None, names=["chr", "start", "end", "junction", "unknown", "strand", "donor_score", "acceptor_score"])
    df = df.sample(n=600)
    
    filtered_score_d = len(df[df['donor_score'] > selected_threshold]['donor_score']) / 600
    filtered_scores_d.append(filtered_score_d)
    filtered_score_a = len(df[df['acceptor_score'] > selected_threshold]['acceptor_score']) / 600
    filtered_scores_a.append(filtered_score_a)
    print(f'{threshold}; filtered_score_d: {filtered_score_d}')
    print(f'{threshold}; filtered_score_a: {filtered_score_a}')

# fig, axes = plt.subplots(1, 2, figsize=(16, 4))  # Adjusting the size for better visibility
fig, axes = plt.subplots(1, 2, figsize=(9, 6))  # Adjusting the size for better visibility

axes[0].set_title(f'Donor Site')
# axes[0].scatter(thresholds_list, filtered_scores_d)
axes[0].plot(thresholds_list, filtered_scores_d, 'o', label='Ratio of Good Donor Site', markersize = 3)

# Calculating the trend line
z = np.polyfit(thresholds_list, filtered_scores_d, 1)  # Fit a 1st degree polynomial to the data
p = np.poly1d(z)  # Create the polynomial equation based on the fit
# Plotting the trend line
axes[0].plot(thresholds_list, p(thresholds_list), "r--", label='Donor trend line', linewidth=3.0)  # "r--" makes the 
axes[0].set_xlabel('Spliced Alignment Count')
axes[0].set_ylabel(f'Ratio of Splice Junctions with Splam Scores above Threshold {selected_threshold}')
axes[0].legend()

axes[1].set_title(f'Acceptor Site')
# axes[1].scatter(thresholds_list, filtered_scores_a)
axes[1].plot(thresholds_list, filtered_scores_a, 'o', label='Ratio of Good Acceptor Site', markersize = 3)

# Calculating the trend line
z = np.polyfit(thresholds_list, filtered_scores_a, 1)  # Fit a 1st degree polynomial to the data
p = np.poly1d(z)  # Create the polynomial equation based on the fit
# Plotting the trend line
axes[1].plot(thresholds_list, p(thresholds_list), "r--", label='Acceptor trend line', linewidth=3.0)  # "r--" makes the 
axes[1].set_xlabel('Spliced Alignment Count')
axes[1].set_ylabel(f'Ratio of Splice Junctions with Splam Scores above Threshold {selected_threshold}')
axes[1].legend()

plt.tight_layout()
plt.subplots_adjust(wspace=0.2)  # Adjust the width of the space between the subplots
plt.savefig('ratio_above_scores.png', dpi=300)