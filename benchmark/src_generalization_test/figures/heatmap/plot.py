import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import itertools
thresholds = [0.95, 0.80, 0.65, 0.50, 0.35, 0.20, 0.05]
db_names = {'GRCm39':'Mouse', 'NHGRI_mPanTro3':'Chimpanzee', 'TAIR10':'Arabidopsis'}
POS_NUM = 12500
NEG_NUM = 25000

def handle_duplicate_names(path):
    # pre: path has '/path/to/file(optversion).ext'
    filename, extension = os.path.splitext(path)
    count = 1
    while os.path.exists(path):
        path = filename + '(' + str(count) + ')' + extension 
        count += 1
    
    return path

def read_inputs(db):
    # positive
    noN_pos_df = pd.read_csv(f'../../pos_test/6_output/combine/averaged_data.noN.{db}.csv')
    noN_pos_df['true_label'] = True
    print(len(noN_pos_df))

    # negative
    noN_neg_df = pd.read_csv(f'../../neg_test/6_output/combine/averaged_data.noN.{db}.csv')
    noN_neg_df['true_label'] = False
    print(len(noN_neg_df))

    noN_pos_df = noN_pos_df.sample(n=POS_NUM, random_state=1091)
    noN_neg_df = noN_neg_df.sample(n=NEG_NUM, random_state=5802)

    noN_merge_df = pd.concat([noN_pos_df, noN_neg_df], axis=0)

    print("noN_pos_df: ", len(noN_pos_df))
    print("noN_neg_df: ", len(noN_neg_df))
    print("noN_merge_df: ", len(noN_merge_df))
    print(noN_merge_df.head())

    return noN_merge_df

def calculate_metrics(df, site, model, threshold):
    true_labels = df['true_label']

    if site == 'donor':
        if model == 'splam':
            scores = df['d_score_splam']
        elif model == 'spliceai':
            scores = df['d_score_spliceai']
    elif site == 'acceptor':
        if model == 'splam':
            scores = df['a_score_splam']
        elif model == 'spliceai':
            scores = df['a_score_spliceai']
    elif site == 'junction': # take minimum between both sites
        if model == 'splam':
            scores = df[['d_score_splam', 'a_score_splam']].min(axis=1)
        elif model == 'spliceai':
            scores = df[['d_score_spliceai', 'a_score_spliceai']].min(axis=1)

    # Apply threshold to scores
    predicted_labels = scores >= threshold
    print(len(predicted_labels), len(true_labels))

    # Get prediction vs. true
    true_positives = np.sum(np.logical_and(predicted_labels, true_labels))
    true_negatives = np.sum(np.logical_and(~predicted_labels, ~true_labels))
    false_positives = np.sum(np.logical_and(predicted_labels, ~true_labels))
    false_negatives = np.sum(np.logical_and(~predicted_labels, true_labels))
    print(true_positives, true_negatives, false_positives, false_negatives)

    # Calculate metrics
    sensitivity = true_positives / (true_positives + false_negatives) # true positive rate
    specificity = true_negatives / (true_negatives + false_positives) # true negative rate
    precision = true_positives / (true_positives + false_positives) # positive predictive value
    npv = true_negatives / (true_negatives + false_negatives) # negative predictive value

    accuracy = (true_positives + true_negatives) / (true_positives + true_negatives + false_positives + false_negatives) # overall correctness

    f1_score = 2 * (precision * sensitivity) / (precision + sensitivity) # harmonic mean of precision and recall

    return sensitivity, specificity, precision, npv, accuracy, f1_score

###### RUNNER ######
def get_stats():

    dbs = ["NHGRI_mPanTro3", "GRCm39", "TAIR10"]
    sites = ['donor', 'acceptor', 'junction']
    models = ['splam', 'spliceai']

    # output
    csv_path = f'./Scores_Data_{POS_NUM}-{NEG_NUM}.csv'
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)

    # generate stats
    full_df = pd.DataFrame(columns=['Database','Site','Model','Sensitivity','Specificity','Precision','NPV','Accuracy','F1_Score', 'Name', 'Threshold'])
    for threshold in reversed(thresholds):
        for db in dbs:
            df = read_inputs(db)
            for site, model in itertools.product(sites, models):
                print(f'Calculating for database: {db}, site: {site}, model: {model}, threshold: {threshold}')
                sensitivity, specificity, precision, npv, accuracy, f1_score = calculate_metrics(df, site, model, threshold)
                full_df.loc[len(full_df)] = [db, site, model, sensitivity, specificity, precision, npv, accuracy, f1_score, db_names[db], threshold]
    
    full_df.reset_index(drop=True, inplace=True)
    print(f'Preview full dataframe\n{full_df}')
    full_df.to_csv(csv_path)

    return full_df

###### PLOTTER ########
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
    fig, axes = plt.subplots(3, 3, figsize=(18, 14), sharex=True, sharey=True)

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

    # Add a common color bar for all subplots
    fig.colorbar(ax.collections[0], ax=axes, location='right', pad=0.05)

    # Save the plot as an image file (PNG format)
    plt.savefig(handle_duplicate_names(f'./Scores_Heatmap_{POS_NUM}-{NEG_NUM}.png'), bbox_inches='tight', dpi=300)

if __name__ == '__main__':
    
    df = get_stats()
    create_plot(df)