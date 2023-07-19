# If we choose a score threshold of 0.8, then the
# (1)sensitivity / (2)accuracy of Splam for (a)donor / (b)acceptor sites / (c)splice junctions  would be
# X%, Y%, and Z% for chimpanzee, mouse, and Arabidopsis

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import itertools
import os
threshold = 0.05

###### INPUT DATA #########
def read_inputs(db):
    # positive
    noN_pos_df = pd.read_csv(f'../../pos_test/6_output/combine/averaged_data.noN.{db}.csv')
    noN_pos_df['true_label'] = True
    print(len(noN_pos_df))

    # negative
    noN_neg_df = pd.read_csv(f'../../neg_test/6_output/combine/averaged_data.noN.{db}.csv')
    noN_neg_df['true_label'] = False
    print(len(noN_neg_df))

    noN_pos_df = noN_pos_df.sample(n=25000, random_state=1091)
    noN_neg_df = noN_neg_df.sample(n=25000, random_state=5802)

    noN_merge_df = pd.concat([noN_pos_df, noN_neg_df], axis=0)

    print("noN_pos_df: ", len(noN_pos_df))
    print("noN_neg_df: ", len(noN_neg_df))
    print("noN_merge_df: ", len(noN_merge_df))
    print(noN_merge_df.head())

    return noN_merge_df

####### CALCULATE METRICS ###########
def calculate_metrics(df, site, model):
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

###### VISUALIZE #######
def plot(df):
    
    # Set the plot style
    sns.set(style="whitegrid")

    # Get unique databases
    databases = df['Database'].unique()

    # Create a plot for each database
    for database in databases:
        # Filter the DataFrame for the current database
        df_database = df[df['Database'] == database]
        
        # Plot using seaborn
        plt.figure(figsize=(8, 5))
        sns.barplot(x='Site', y='F1_Score', hue='Model', data=df_database, palette='Set2')
        
        # Set plot title and labels
        plt.title(f'F1 Scores for {database}')
        plt.xlabel('Site')
        plt.ylabel('F1 Score')
        
        # Show the legend
        plt.legend(title='Model')
        plt.yscale('log')
        
        # Save the plot as an image file (PNG format)
        plt.savefig(f'./{threshold:.1e}/F1_Scores_{database}.png', bbox_inches='tight', dpi=300)
    


###### RUNNER ######
def run():

    dbs = ["NHGRI_mPanTro3", "GRCm39", "TAIR10"]
    sites = ['donor', 'acceptor', 'junction']
    models = ['splam', 'spliceai']

    # output
    csv_path = f'./{threshold:.1e}/result.csv'
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)

    # generate and print statistics
    with open(csv_path, 'w') as f:
        f.write('Database,Site,Model,Sensitivity,Specificity,Precision,NPV,Accuracy,F1_Score\n')
        for db in dbs:
            df = read_inputs(db)
            for site, model in itertools.product(sites, models):
                print('-'*120)
                print(f'Calculating for database: {db}, site: {site}, model: {model}')
                sensitivity, specificity, precision, npv, accuracy, f1_score = calculate_metrics(df, site, model)
                print(f'\tSensitivity: {sensitivity}\n\tSpecificity: {specificity}\n\tPrecision: {precision}\n\tNPV: {npv}\n\tAccuracy: {accuracy}\n\tF1 Score: {f1_score}\n')

                f.write(f'{db},{site},{model},{sensitivity},{specificity},{precision},{npv},{accuracy},{f1_score}\n')
    
    # create plot
    print('Plotting figure...')
    stats = pd.read_csv(csv_path)
    plot(stats)

if __name__ == '__main__':
    run()
