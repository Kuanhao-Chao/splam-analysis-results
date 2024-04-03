import os
import numpy as np
import matplotlib.pyplot as plt

project_root = "/ccb/cybertron/khchao/splam-analysis-results/"
# albation_study_root=f'{project_root}results/albation_study/MODEL/subset_10000/'
albation_study_root=f'{project_root}train/src/MODEL/'
output_dir = f'{albation_study_root}img/'
os.makedirs(output_dir, exist_ok=True)

res_type = ['TRAIN', 'TEST']

def moving_average(data, window_size):
    """Compute the moving average of a list Xof numbers."""
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')

def calculate_f1(precision, recall):
    """Calculate F1 score from precision and recall."""
    return [2 * p * r / (p + r) if (p + r) else 0 for p, r in zip(precision, recall)]

for eval_target in ['A', 'D', 'J']:
    if eval_target == 'A':
        metrics = [
            "_A_accuracy.txt",
            "_A_auprc.txt",
            "_A_threshold_precision.txt",
            "_A_threshold_recall.txt",
        ]
    elif eval_target == 'D':
        metrics = [
            "_D_accuracy.txt",
            "_D_auprc.txt",
            "_D_threshold_precision.txt",
            "_D_threshold_recall.txt",
        ]
    elif eval_target == 'J':
        metrics = [
            "_J_threshold_f1.txt",
            "_J_threshold_precision.txt",
            "_J_threshold_recall.txt",
        ]

    for window_size in [1, 5, 20, 40, 80]:
        # Initialize a dictionary to hold all scores
        scores_dict = {}
        for res in res_type:
            for metric in metrics:
                scores_dict[f'{res}_{metric}'] = {}
                print(f'{res}; {metric}')
                for seq_len in [40, 100, 200, 400, 600, 800]:
                    key = f'Splam_{seq_len}'
                    # scores_key = f'{exp}_{rs_value}'
                    key_dir = f'{albation_study_root}{key}/'
                    file_path = f'{key_dir}LOG/{res}/{res.lower()}{metric}'
                    print("file_path: ", file_path)
                    if os.path.exists(file_path):
                        with open(file_path, 'r') as file:
                            lines = file.readlines()
                            scores = [float(line.strip()) for line in lines]
                            # Initialize a dictionary if it does not exist
                            if key not in scores_dict[f'{res}_{metric}']:
                                scores_dict[f'{res}_{metric}'][key] = []
                            scores_dict[f'{res}_{metric}'][key].append(scores)
                    else:
                        print(f'{file_path} does not exist.')

        if eval_target == 'J':
            # After collecting precision and recall, calculate F1 for J
            for res in res_type:  # Ensure `res` is still accessible
                precision_key = f'{res}__J_threshold_precision.txt'
                recall_key = f'{res}__J_threshold_recall.txt'
                f1_key = f'{res}__J_threshold_f1.txt'
                for scores_key in scores_dict[precision_key].keys():  # Iterate through each experiment condition
                    precision = scores_dict[precision_key][scores_key]
                    recall = scores_dict[recall_key][scores_key]
                    f1_scores = [calculate_f1(p, r) for p, r in zip(precision, recall)]
                    scores_dict[f1_key] = scores_dict.get(f1_key, {})
                    scores_dict[f1_key][scores_key] = f1_scores

        # Visualization
        fig, axs = plt.subplots(len(res_type), len(metrics), figsize=(40,12), dpi=800)  # Adjusted figsize and added dpi for clarity
        for i, res in enumerate(res_type):
            for j, metric in enumerate(metrics):
                ax = axs[i, j]
                data_dict = scores_dict[f'{res}_{metric}']
                if data_dict:
                    for scores_key, data in data_dict.items():
                        print("scores_key: ", scores_key)
                        if len(data) > 0:
                            for scores in data:
                                if len(scores) > window_size:
                                    smoothed_scores = moving_average(scores, window_size)
                                    ax.plot(smoothed_scores, label=scores_key)
                                else:
                                    ax.plot(scores, label=scores_key)  # Plot raw data if too short for smoothing
                    ax.set_title(f'{res} {metric}', fontsize=27)  # Larger title
                    ax.set_xlabel('X-axis Label', fontsize=16)  # X-axis label size
                    ax.set_ylabel('Y-axis Label', fontsize=16)  # Y-axis label size
                    ax.tick_params(axis='both', which='major', labelsize=12)  # Adjust tick size for both axes
                    ax.legend(fontsize=12)  # Adjust legend font size
                    # ax.set_ylim([0, 1])  # Set y-axis scale from 0 to 1
                else:
                    ax.set_title(f'{res}_{metric} (No Data)', fontsize=16)  # Larger title
                    ax.axis('off')
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.4, wspace=0.2)  # Add more space between subplots
        plt.savefig(f'{output_dir}seqlen_vary_smooth{window_size}_{eval_target}.png', dpi=100)