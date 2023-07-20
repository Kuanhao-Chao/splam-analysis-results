# run in /benchmark/src_compare_results/

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os

def run_plotter(type):

    # obtain the averaged data from POSITIVE sample
    for db in ['GRCm39', 'Mmul_10', 'NHGRI_mPanTro3', 'TAIR10']:
        avg_df = pd.read_csv(f'../../pos_test/6_output/combine/averaged_data.{type}.{db}.csv')
        avg_df = revive_list(avg_df, ['spliceai_version'], ['int'])
        get_stats(avg_df, type, db)

    avg_df = pd.read_csv(f'../../pos_test/6_output/combine/averaged_data.{type}.TAIR10.csv')
    avg_df = revive_list(avg_df, ['spliceai_version'], ['int'])
    tair(avg_df, type, db)

def revive_list(df, col_names, dtypes):
    
    for name, dt in zip(col_names, dtypes):
        # revive the lists in column from string representations into np arrays
        df[name] = df[name].apply(lambda x: np.array(x.strip('[]').split(','), dtype=dt))

    return df

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
    plt.savefig(fig_path, dpi=300)
    print(f'Saved figure to {fig_path}.')

def get_stats(avg_df, type, db):

    # getting chrom sizes
    chrs = {}
    with open(f'../../data/{db}_assembly_report.txt', 'r') as file:       
        # read the file line by line
        for line in file:  
            if line.startswith('#'):
                continue
            # split by tabs
            columns = line.strip().split('\t')
            refseq_name = columns[6]
            chrom_size = int(columns[8])

            # store the key-value pair in the dictionary
            chrs[refseq_name] = chrom_size
    
    # get mean statistics 
    df = avg_df.drop(['name', 'expected_score', 'strand', 'd_dimer', 'a_dimer', 'spliceai_version'], axis=1)
    df['ACTG Content'] = df.iloc[:,7:11].sum(axis=1)
    df['N Content'] = df.iloc[:,11:18].sum(axis=1)
    df.drop(df.columns[7:18], axis=1, inplace=True)
    df['Intron Length'] = df['end'] - df['start']
    df['dist_to_start'] = df['start'] - 0
    df['dist_to_end'] = df['seqid'].map(chrs) - df['end']
    df['Distance to Ending'] = df[['dist_to_start', 'dist_to_end']].min(axis=1)
    df['ACTG Content'] = df['ACTG Content'] / (df['Intron Length'] + 10400)
    df['N Content'] = df['N Content'] / (df['Intron Length'] + 10400)
    mean_stats_df = df[['Intron Length', 'Distance to Ending', 'ACTG Content', 'N Content']].mean(axis=0).astype('float64')

    mean_stats_df.to_csv(f'mean_stats.{type}.{db}.csv')

    print(df, '\nStats:\n', mean_stats_df)


'''Investigating TAIR scores in SpliceAI vs. Splam'''
def tair(avg_df, type, db):

    # getting chrom sizes
    chrs = {}
    with open(f'../../data/{db}_assembly_report.txt', 'r') as file:       
        # read the file line by line
        for line in file:  
            if line.startswith('#'):
                continue
            # split by tabs
            columns = line.strip().split('\t')
            refseq_name = columns[6]
            chrom_size = int(columns[8])

            # store the key-value pair in the dictionary
            chrs[refseq_name] = chrom_size

    # base statistics
    df = avg_df.drop(['name', 'expected_score', 'strand', 'd_dimer', 'a_dimer', 'spliceai_version'], axis=1)
    df['ACTG Content'] = df.iloc[:,7:11].sum(axis=1)
    df['N Content'] = df.iloc[:,11:18].sum(axis=1)
    df.drop(df.columns[7:18], axis=1, inplace=True)
    df['Intron Length'] = df['end'] - df['start']
    df['dist_to_start'] = df['start'] - 0
    df['dist_to_end'] = df['seqid'].map(chrs) - df['end']
    df['Distance to Ending'] = df[['dist_to_start', 'dist_to_end']].min(axis=1)
    df['ACTG Content'] = df['ACTG Content'] / (df['Intron Length'] + 10400)
    df['N Content'] = df['N Content'] / (df['Intron Length'] + 10400)


    d_df = pd.melt(df, value_vars=['d_score_spliceai', 'd_score_splam'], var_name='Method', value_name='Score', 
                   id_vars=['N Content', 'Intron Length', 'Distance to Ending']).replace('d_score_spliceai', 'SpliceAI').replace('d_score_splam', 'SPLAM')
    a_df = pd.melt(df, value_vars=['a_score_spliceai', 'a_score_splam'], var_name='Method', value_name='Score', 
                   id_vars=['N Content', 'Intron Length', 'Distance to Ending']).replace('a_score_spliceai', 'SpliceAI').replace('a_score_splam', 'SPLAM')
    print('Donor and Acceptor Stats\n\n', d_df, '\n', a_df)

    # investigating intron len specifically lower end
    threshold = 400
    for above in [True, False]:
        for donor in [True, False]:
            if donor:
                df = d_df
                site = 'Donor'
            else:
                df = a_df
                site = 'Acceptor'
            if above:
                df = df[df['Intron Length'] > threshold]
                disc = '>'
            else:
                df = df[df['Intron Length'] < threshold]
                disc = '<'

            sns.set(font_scale=0.8)
            sns.set_palette('deep')
            plot_params = {'s': 8, 'alpha': 0.6}
            sns.jointplot(data=df, x='Intron Length', y='Score', hue='Method', height=8, ratio=4, **plot_params)
            plt.title(f'{site} Intron Length ({disc}400) vs. Score for {db} Dataset')
            save_fig(f'./{site}{disc}{threshold}_jointplot_len_vs_score.{type}.{db}.png')


if __name__ == '__main__':

    type = 'noN'
    run_plotter(type)