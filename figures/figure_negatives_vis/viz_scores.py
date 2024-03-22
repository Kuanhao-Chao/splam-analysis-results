import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from io import StringIO

# for log in [False, True]:
for log in [False]:
    for thresholds in range(2, 101):
        neg_pred_file = f'/ccb/cybertron/khchao/splam-analysis-results/src_tools_evaluation/splam_result/RELEASED/neg_{thresholds}/LOG/_junction_score.bed'
        if log:
            output_png = f'/ccb/cybertron/khchao/splam-analysis-results/src_tools_evaluation/splam_result/RELEASED/neg_{thresholds}/LOG/_junction_score_log.png'
        else:
            output_png = f'/ccb/cybertron/khchao/splam-analysis-results/src_tools_evaluation/splam_result/RELEASED/neg_{thresholds}/LOG/_junction_score.png'

        df = pd.read_csv(neg_pred_file, sep="\t", header=None, names=["chr", "start", "end", "junction", "unknown", "strand", "donor_score", "acceptor_score"])
        # Creating the subplots: 1 row, 2 columns
        fig, axes = plt.subplots(1, 2, figsize=(12, 3))  # Adjusting the size for better visibility

        # Plotting the donor score distribution on the first subplot
        sns.histplot(df['donor_score'], ax=axes[0], color="blue", label='Donor Score', kde=True, bins=100)
        if log:
            axes[0].set_yscale('log')
        axes[0].set_title('Distribution of Donor Scores')
        axes[0].set_xlabel('Score')
        axes[0].set_ylabel('Count')
        axes[0].legend()

        # Plotting the acceptor score distribution on the second subplot
        sns.histplot(df['acceptor_score'], ax=axes[1], color="red", label='Acceptor Score', kde=True, bins=100)
        if log:
            axes[1].set_yscale('log')
        axes[1].set_title('Distribution of Acceptor Scores')
        axes[1].set_xlabel('Score')
        axes[1].set_ylabel('Count')
        axes[1].legend()

        plt.tight_layout()
        plt.savefig(output_png, dpi=300)
        # plt.show()
