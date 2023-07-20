import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os

def handle_duplicate_names(path):
    # pre: path has '/path/to/file(optversion).ext'
    filename, extension = os.path.splitext(path)
    count = 1
    while os.path.exists(path):
        path = filename + '(' + str(count) + ')' + extension 
        count += 1
    
    return path

def save_fig(figpath):
    fig_path = handle_duplicate_names(figpath)
    os.makedirs(os.path.dirname(fig_path), exist_ok=True)
    plt.savefig(fig_path, dpi=500)
    print(f'Saved figure to {fig_path}.')

def plot_hist(type, dataset):
    # defining paper figure for comparing SPLAM and SpliceAI in histogram across all species for donor and acceptor
    
    # data gathering from all species
    species = ['GRCm39', 'Mmul_10', 'NHGRI_mPanTro3', 'TAIR10']
    sites = ['d_splam', 'd_spliceai', 'a_splam', 'a_spliceai']
    cols = [x+'_'+y for x in species for y in sites]
    all_species_df = pd.DataFrame()

    for idx, db in zip(range(0, len(cols), 4), species):
        avg_ofp = f'../../{dataset}_test/6_output/combine/averaged_data.{type}.{db}.csv'
        df = pd.read_csv(avg_ofp)

        t_df = pd.DataFrame(columns=cols[idx:idx+4])
        t_df.iloc[:,0] = df['d_score_splam']
        t_df.iloc[:,1] = df['d_score_spliceai']
        t_df.iloc[:,2] = df['a_score_splam']
        t_df.iloc[:,3] = df['a_score_spliceai']
        
        all_species_df = pd.concat([all_species_df, t_df], axis=1)
    print(all_species_df)

    # plot setup
    sns.set(font_scale=0.8)
    sns.set_palette('deep')
    f, axs = plt.subplots(2, 4, figsize=(16, 10), sharey=True, sharex=True)
    plt.subplots_adjust(hspace=None,wspace=None)
    f.supylabel('Density', size=17, weight=400)
    f.supxlabel('Score', size=17, weight=400)

    # calculate quartiles
    quartiles = all_species_df.quantile([0.25, 0.50, 0.75], axis=0)
    diff = pd.DataFrame(quartiles.loc[0.75] - quartiles.loc[0.25]).T
    quartiles = pd.concat([quartiles, diff])
    print(quartiles)
    quartiles.to_csv(f'./{dataset}_quartiles_{type}.csv')

    # traverse subplots column-major (top-down, left-right)
    axs = axs.flatten(order='F')
    for i in range(0,16,2):

        # plot first to get lines
        data = all_species_df.iloc[:,i:i+2]
        ax = axs[i//2]
        lines = sns.kdeplot(data=data, ax=ax, clip=(0.0, 1.0)).get_lines()
        x1, y1 = lines[0].get_data() # SpliceAI - orange
        x2, y2 = lines[1].get_data() # SPLAM - blue

        # plot again with fill and format
        p = sns.kdeplot({'splam': data.iloc[:,0], 'SpliceAI': data.iloc[:,1]}, ax=ax, clip=(0.0, 1.0), fill=True, alpha=0.35, label=["SpliceAI", "splam"])
        sns.move_legend(p, 'upper left')
        ax.tick_params(axis='x', labelsize=10)

        # plot quartiles
        quarts = quartiles.iloc[:,i:i+2].values.flatten()
        print(quarts)
        for j in range(0,6,2): 
            splam_max = np.interp(quarts[j], x2, y2)
            spliceai_max = np.interp(quarts[j+1], x1, y1)
            ax.vlines(quarts[j], 0, splam_max, colors='b', alpha=0.4, linestyles='dashed')
            ax.vlines(quarts[j+1], 0, spliceai_max, colors='r', alpha=0.4, linestyles='dashed')
            if j == 2: # median
                ax.vlines(quarts[j], 0, splam_max, colors='navy', alpha=1, linestyles='solid')
                ax.vlines(quarts[j+1], 0, spliceai_max, colors='r', alpha=1, linestyles='solid')
    
        if i % 2 == 0 and i < 8:
            # x-axis species label
            axs[i].set_title(species[i//2], size=14, weight=500)

    # y-axis site label
    axs[0].set_ylabel('Donor Site', size=14, weight=500)
    axs[1].set_ylabel('Acceptor Site', size=14, weight=500)
    axs[0].tick_params(axis='y', labelsize=12, width=12)
    axs[1].tick_params(axis='y', labelsize=12, width=12)

    # display and save
    plt.tight_layout(rect=(0, 0, 1, 0.95))
    save_fig(f'./{dataset}_scores_compare.{type}.png')

if __name__ == '__main__':

    type = 'noN'
    for dataset in ['pos', 'neg']:
        plot_hist(type, dataset)