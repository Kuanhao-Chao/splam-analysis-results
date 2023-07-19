
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import itertools

def read_input():

    # consts
    thresholds = [0.1, 0.01, 0.001, 0.00001, 0.00001, 0.000001]
    db_names = {'GRCm39':'Mouse', 'NHGRI_mPanTro3':'Chimpanzee', 'TAIR10':'Arabidopsis'}

    # read in the files 
    full_df = pd.DataFrame(columns=['Database','Site','Model','Sensitivity','Specificity','Precision','NPV','Accuracy','F1_Score','Name','Threshold'])
    for threshold in reversed(thresholds):
        df = pd.read_csv(f'../f1_scores/{threshold:.0e}/result.csv')
        df['Name'] = df['Database'].map(db_names)
        df['Threshold'] = f'{threshold}'
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
    df.drop(columns=['Model','Site'], inplace=True)
    print(score_df)

    # create the plots

    pass

if __name__ == '__main__':

    df = read_input()
    create_plot(df)
    
