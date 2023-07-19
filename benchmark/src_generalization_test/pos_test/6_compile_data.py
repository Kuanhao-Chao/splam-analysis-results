import os
from progress.bar import Bar
import numpy as np
import pandas as pd
from pyfaidx import Fasta 
import itertools

'''Run all parts in sequence'''
def run_all(db_name, overwrite_all=False):
    types = ['noN', 'N']
    versions = [1, 2, 3, 4, 5]

    print('AGGREGATING SPLAM DATA...')
    splam_aggregator(db_name, overwrite_all)
    print('SPLAM DONE.')

    print('AGGREGATING SPLICEAI DATA...')
    for type, version in itertools.product(types, versions):
        spliceai_aggregator(type, version, db_name, overwrite_all)
    print('SPLICEAI DONE.')

    print('COMBINING...')
    for type in types:
        combine(type, db_name, overwrite_all)
    print('COMBINE DONE.')


#######################
# RUN SPLAM PORTION
#######################
def splam_aggregator(db_name, overwrite=False):

    # collect the data to be aggregated from multiple sources
    splam_score_file = f'./4_output/{db_name}/score.bed'
    output_file = f'./6_output/Splam/{db_name}.splam_data.csv'

    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # make csv file if it does not exist, write to pd dataframe
    if not os.path.exists(output_file) or overwrite:
        collect_splam(splam_score_file, db_name, output_file)

'''Collect the Splam scores into a single dataframe'''
def collect_splam(score_file, db_name, output_file):
    # read score data from score file
    df = pd.read_csv(score_file, sep='\t', header=None, \
        names=['chrom', 'chromStart(donor)', 'chromEnd(acceptor)', 'name', 'score', 'strand', 'donorScore', 'acceptorScore'])

    # add new columns for the dimers
    df['donorDimer'] = ''
    df['acceptorDimer'] = ''

    # ping the fasta file using pyfaidx to obtain the sequences as a dictionary
    genes = Fasta(f'../data/{db_name}_genomic.fa', sequence_always_upper=True)
    print(f'Found {len(genes.keys())} unique genes in {db_name} fasta file.')

    # iterate over every coordinate and extract corresponding dimer from fasta
    pbar = Bar('[Info] Getting dimers...', max=len(df))
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    for index, row in df.iterrows():
        # obtain the start and end sequences from the fasta file
        chromosome = row['chrom']
        strand = row['strand']
        donor_start = int(row['chromStart(donor)'])
        acceptor_end = int(row['chromEnd(acceptor)'])

        # extract the dimers from the fasta file
        if strand == '+': 
            donor_dimer = str(genes[chromosome][donor_start:donor_start+2])
            acceptor_dimer = str(genes[chromosome][acceptor_end-2:acceptor_end])

        elif strand == '-':
            acceptor_dimer = str(genes[chromosome][donor_start:donor_start+2].reverse.complement)
            donor_dimer = str(genes[chromosome][acceptor_end-2:acceptor_end].reverse.complement)

        # insert the dimers into the pandas dataframe
        df.at[index, 'donorDimer'] = donor_dimer
        df.at[index, 'acceptorDimer'] = acceptor_dimer

        pbar.next()

    pbar.finish()
    print(df.head())

    # convert df to csv file 
    df.to_csv(output_file, index=False)
    print(f'Full data csv file saved to {output_file}')


#######################
# RUN SPLICEAI PORTION
#######################
def spliceai_aggregator(TYPE, SPLICEAI_VERSION, db, overwrite=False):

    # define identifiers for this run
    print('*'*140)
    print(f'Parsing for type {TYPE}, SpliceAI version {SPLICEAI_VERSION}, database {db}')

    # input filepaths
    dirpath = f'./5_output/spliceai_result_{SPLICEAI_VERSION}/{db}/'
    d_score_file = f'{dirpath}spliceai_all_seq.score.d.{TYPE}.{db}.tsv'
    a_score_file = f'{dirpath}spliceai_all_seq.score.a.{TYPE}.{db}.tsv'
    name_file = f'{dirpath}spliceai_all_seq.name.{TYPE}.{db}.tsv'
    splam_data = f'./6_output/Splam/{db}.splam_data.csv'

    # output filepaths
    id = f'.v{SPLICEAI_VERSION}.{TYPE}.{db}'
    csv_path = f'./6_output/SpliceAI/{db}/spliceai_data{id}.csv'
    donor_path = f'./6_output/SpliceAI/{db}/donor_data{id}.tsv'
    acceptor_path = f'./6_output/SpliceAI/{db}/acceptor_data{id}.tsv'
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    
    # collect the full data into a dataframe
    if not os.path.exists(csv_path) or overwrite:
        score_df = collect_spliceai(d_score_file, a_score_file, donor_path, acceptor_path)
        full_df = index_compare(splam_data, name_file, score_df)
        write_df(full_df, csv_path)

'''Parse the SpliceAI donor and acceptor files into a single dataframe'''
def collect_spliceai(d_file, a_file, donor_path, acceptor_path):

    # read donor and acceptor spliceai files into list
    print('Reading donor score file...')
    with open(d_file, 'r') as d_read:
        d_lines = d_read.read().splitlines()
    print('Reading acceptor score file...')
    with open(a_file, 'r') as a_read:
        a_lines = a_read.read().splitlines()

    # collect 400 scores on donor side, then donor position score for each intron
    pbar = Bar('Collecting donor scores...', max=len(d_lines))
    donors = []
    donor_scores = []
    for line in d_lines:
        values = line.strip().split()
        donor_row = values[:400]
        donors.append(donor_row)
        donor_scores.append(donor_row[199])
        pbar.next()
    pbar.finish()

    # collect 400 scores on acceptor side, then acceptor position score for each intron
    pbar = Bar('Collecting acceptor scores...', max=len(a_lines))
    acceptors = []
    acceptor_scores = []
    for line in a_lines:
        values = line.split()
        acceptor_row = values[-400:]
        acceptors.append(acceptor_row)
        acceptor_scores.append(acceptor_row[200])
        pbar.next()
    pbar.finish()

    # compile scores into a single dataframe
    score_df = pd.DataFrame()
    score_df['donor_score'] = pd.Series(donor_scores, dtype='float64')
    score_df['acceptor_score'] = pd.Series(acceptor_scores, dtype='float64')
    donor_df = pd.DataFrame(donors, dtype='float64')
    acceptor_df = pd.DataFrame(acceptors, dtype='float64')

    # sanity check 
    assert(len(score_df) == 25000)
    assert(len(donor_df) == 25000)
    assert(len(acceptor_df) == 25000)

    # save the donor and acceptor scores first
    #donor_df.to_csv(donor_path, sep='\t', header=None, index=0) 
    #acceptor_df.to_csv(acceptor_path, sep='\t', header=None, index=0) 

    # print results of score 
    print(f'Preview SpliceAI scores:\n{score_df}')
    print(f'Averages:\n{score_df.mean()}')
    print('Finished parsing')

    return score_df

'''Compare the indices between the SPLAM data file and SpliceAI output and get only the matches'''
def index_compare(full_data_file, name_file, score_df):

    # read files and create dfs (df1 = SPLAM file; df2 = SpliceAI file w. scores)
    splam_full_df = pd.read_csv(full_data_file)
    bed_names = ['chromosome', 'start', 'end', 'strand']
    bed_types = {'chromosome': str, 'start': int, 'end': int, 'strand': str} 
    spliceai_name_df = pd.read_csv(name_file, delimiter=' ', header=None, names=bed_names, dtype=bed_types)

    # identify the dataframes to inner merge (df1 and df2)
    df1 = splam_full_df
    df2 = pd.concat([spliceai_name_df, score_df], axis=1)
    df2.columns = ['id', 'start', 'end', 'strand', 'd_score_spliceai', 'a_score_spliceai']

    # get the identifiers (seqid, start, end) as a list of tuples -> show uniques for each 
    print('Filtering full data...')
    df1_rows = df1.iloc[:,:3].apply(tuple,axis=1)
    df2_rows = df2.iloc[:,:3].apply(tuple,axis=1)
    print(f'\t\tRows\tUniques\nSPLAM:\t\t{len(df1_rows)}\t{len(df1_rows.drop_duplicates())}\nSpliceAI:\t{len(df2_rows)}\t{len(df2_rows.drop_duplicates())}')
    df1.drop_duplicates(subset=['chrom', 'chromStart(donor)', 'chromEnd(acceptor)', 'strand'], inplace=True)

    # merge the two dataframes by ID, sort result, and drop duplicates 
    # NOTE: this will still keep some duplicate results as SpliceAI may give different scores to same sequence 
    filtered_df = df1.merge(df2, how='inner', left_on=['chrom', 'chromStart(donor)', 'chromEnd(acceptor)', 'strand'], 
                            right_on=['id', 'start', 'end', 'strand'], sort=True)
    
    assert(len(filtered_df) == 25000)

    print(f'Merged rows: {len(filtered_df)}')
    print(f'Preview filtered SPLAM values:\n{filtered_df}')
    print('Finished filtering')

    return filtered_df

'''Combine the full data with SpliceAI scores and save to csv file'''
def write_df(full_df, csvpath):

    # drop unneccesary columns and fix indices
    full_df.drop(['id', 'start', 'end'], inplace=True, axis=1)
    full_df.reset_index(drop=True, inplace=True)

    # reorder and rename the final df columns
    full_df = full_df.reindex(columns=['chrom', 'chromStart(donor)', 'chromEnd(acceptor)', 'name', 'score', 'strand', 'donorDimer', 'acceptorDimer', 
                                       'donorScore', 'acceptorScore', 'd_score_spliceai', 'a_score_spliceai'])
    full_df.columns = ['seqid', 'start', 'end', 'name', 'expected_score', 'strand', 'd_dimer', 'a_dimer', 
                       'd_score_splam', 'a_score_splam', 'd_score_spliceai', 'a_score_spliceai']

    print(f'Preview final saved df:\n{full_df}')

    # save to file
    os.makedirs(os.path.dirname(csvpath), exist_ok=True)
    full_df.to_csv(csvpath, index=False)
    print(f'Saved csv file to {csvpath}.')


#######################
# COMBINE BOTH PARTS
#######################
def combine(type, db, overwrite=False):

    # identifier
    id = [type, db]
    print(f'Running for type: {type}, database: {db}')

    # obtain the aggregated and averaged data from all 5 versions of spliceai and write to file
    agg_ofp = f'./6_output/combine/aggregate_data.{type}.{db}.csv'
    avg_ofp = f'./6_output/combine/averaged_data.{type}.{db}.csv'
    if not os.path.exists(avg_ofp) or overwrite:
        merge_df = aggregate_data(id)
        agg_df = get_counts(merge_df, id, agg_ofp)
        get_averages(agg_df, avg_ofp)

'''Gets the aggregated data from all 5 SpliceAI versions'''
def aggregate_data(id):

    # define the df which will collect all scores together
    aggregate_df = pd.DataFrame(columns=['seqid', 'start', 'end', 'name', 'expected_score', 'strand', 'd_dimer', 'a_dimer', 
                                            'd_score_splam', 'a_score_splam', 'd_score_spliceai', 'a_score_spliceai', 'spliceai_version'])

    # aggregate all versions of spliceai
    for version_num in range(1, 6):
    
        # read the input file
        ifp = f'./6_output/SpliceAI/{id[1]}/spliceai_data.v{version_num}.{id[0]}.{id[1]}.csv'
        full_df = pd.read_csv(ifp)

        # sanity check
        assert(len(full_df) == 25000)

        # add column indicating which version of the model data is from
        full_df['spliceai_version'] = version_num

        # concat onto aggregate df
        aggregate_df = pd.concat([aggregate_df, full_df], axis=0)

    # group the dataframe by common columns and aggregate the spliceai scores into a list
    merged_df = aggregate_df.groupby(['seqid', 'start', 'end', 'name', 'expected_score', 'strand', 'd_dimer', 'a_dimer', 
                                      'd_score_splam', 'a_score_splam']).agg(lambda x: x.tolist()).reset_index()

    print(f'Preview merged_df:\n{merged_df}')
    
    # sanity check
    assert(len(merged_df) == 25000)
    for i, row in merged_df.iterrows():
        assert(len(row['d_score_spliceai']) == 5)
        assert(len(row['a_score_spliceai']) == 5)

    return merged_df

'''Gets the counts of the individual nucleotides in the input sequences used by SpliceAI'''
def get_counts(df, id, ofp):

    # get the fasta file with the sequence used by spliceai
    fa_file = f'./3_output/{id[1]}/seq_{id[0]}.fa'
    fa = Fasta(fa_file, key_function = lambda x: x[:-2], duplicate_action='first')
    
    # iterate over each row and get the nucleotide counts
    l_set = ['A', 'C', 'T', 'G', 'N', 'R', 'Y', 'K', 'M', 'S', 'W']
    counts = pd.DataFrame(columns=l_set)
    pbar = Bar('Getting counts...', max=len(df))
    for i, row in df.iterrows():

        # the key used to search for the relevant sequence
        idx = f"{row['seqid']};{row['start']};{row['end']}"

        # this is the sequence from start to end, with 5.2k flanking seq on each side
        sequence = fa[idx][:].seq
        length = row['end'] - row['start'] + 10400
        assert(len(sequence) == length)

        # get the counts
        l = [letter for letter in sequence]
        d = {k:l.count(k) for k in l_set}
        count = pd.DataFrame([d])
        counts = pd.concat([counts, count], axis=0)

        pbar.next()
    pbar.finish()

    counts.reset_index(drop=True, inplace=True)
    df = pd.concat([df, counts], axis=1)

    print(f'Preview of df with counts:\n{df}', flush=True)

    # save to csv
    os.makedirs(os.path.dirname(ofp), exist_ok=True)
    df.to_csv(ofp, index=False)
    print(f'Saved agg with counts csv file to {ofp}.')

    return df

'''Applies averages across all scores'''
def get_averages(df, ofp):

    # take the averages of the spliceai score columns
    df['d_score_spliceai'] = df['d_score_spliceai'].apply(np.average)
    df['a_score_spliceai'] = df['a_score_spliceai'].apply(np.average)

    print(f'Preview df with averages:\n{df}')
    
    # save to csv
    os.makedirs(os.path.dirname(ofp), exist_ok=True)
    df.to_csv(ofp, index=False)
    print(f'Saved average csv file to {ofp}.')

    return df


if __name__ == '__main__':

    if os.getcwd() != 'pos_test':
        os.chdir('/home/smao10/splam-analysis-results/benchmark/src_generalization_test/pos_test/')

    dbs = ['GRCm39', 'Mmul_10', 'NHGRI_mPanTro3', 'TAIR10']
    nums = [0,1,2,3] # CHANGEME

    for num in nums:
        run_all(dbs[num])