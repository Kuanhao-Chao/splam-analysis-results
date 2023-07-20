
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import itertools
thresholds = [0.95, 0.80, 0.65, 0.50, 0.35, 0.20, 0.05]
db_names = {'GRCm39':'Mouse', 'NHGRI_mPanTro3':'Chimpanzee', 'TAIR10':'Arabidopsis'}

def handle_duplicate_names(path):
    # pre: path has '/path/to/file(optversion).ext'
    filename, extension = os.path.splitext(path)
    count = 1
    while os.path.exists(path):
        path = filename + '(' + str(count) + ')' + extension 
        count += 1
    
    return path

def read_input():

    # read in the files 
    full_df = pd.DataFrame(columns=['Database','Site','Model','Sensitivity','Specificity','Precision','NPV','Accuracy','F1_Score','Name','Threshold'])
    for threshold in reversed(thresholds):
        df = pd.read_csv(f'../f1_scores/{threshold:.1e}/result.csv')
        df['Name'] = df['Database'].map(db_names)
        df['Threshold'] = threshold
        full_df = pd.concat([full_df, df], axis=0)
        
    full_df.reset_index(drop=True, inplace=True)

    print(f'Preview full dataframe\n{full_df}')

    return full_df

def create_plot(df):

    # process the input
    df.drop(['Specificity','NPV','F1_Score'], axis=1, inplace=True)
    score_df = pd.melt(df, id_vars=['Database','Site','Model','Name','Threshold'], var_name='Score_Type', value_name='Score')
    combine_model_site = lambda row : f"{row['Model']}-{row['Site']}"
    score_df['Model-Site'] = score_df.apply(combine_model_site, axis=1)
    score_df.drop(columns=['Database','Model','Site'], inplace=True)
    df = score_df.reindex(columns=['Name','Score_Type','Model-Site','Threshold','Score'])
    print(df)

    # Get unique values for 'Model-Site' and 'Threshold'
    names = df['Name'].unique()
    score_types = df['Score_Type'].unique()

    # Create a 3x3 grid of subplots
    fig, axes = plt.subplots(3, 3, figsize=(18, 14), sharex=False, sharey=True)

    # Loop through each subplot and plot the heatmap
    for i, score_type in enumerate(score_types):
        for j, name in enumerate(names):
            ax = axes[i, j]
            subset = df[(df['Name'] == name) & (df['Score_Type'] == score_type)].reset_index(drop=True)
            subset['Score'] *= 100
            heatmap_subset = subset.pivot(columns='Threshold', index='Model-Site', values='Score')
            print(heatmap_subset)
            sns.heatmap(heatmap_subset, ax=ax, annot=True, cmap=sns.color_palette('flare', as_cmap=True), cbar=False, fmt='3.1f', vmin=10.0, vmax=100.0, xticklabels=['{:.2f}'.format(i) for i in reversed(thresholds)])
            #ax.set_title(f'{name} - {score_type}')
            ax.set_xlabel('')
            ax.set_ylabel('')

    # Set labels for the x and y axes of the last row and last column of subplots
    # for ax, name in zip(axes[0], names):
    #     ax.set_xlabel(f'{name}', fontsize=15, labelpad=3)
    # for ax, type in zip(axes[:, 0], score_types):
    #     ax.set_ylabel(f'{type}', fontsize=15, labelpad=3)
    
    # fig.supxlabel('Threshold       ', fontsize=18)
    # fig.supylabel('Model & Site', fontsize=18)
    # fig.suptitle('Collected Scores for Different Species, Models, and Sites at Various Thresholds         ', fontsize=21)

    # Add a common color bar for all subplots
    fig.colorbar(ax.collections[0], ax=axes, location='right', pad=0.05)

    # Save the plot as an image file (PNG format)
    plt.savefig(handle_duplicate_names(f'./Scores_Heatmap.png'), bbox_inches='tight', dpi=300)

if __name__ == '__main__':
    
    df = read_input()
    create_plot(df)