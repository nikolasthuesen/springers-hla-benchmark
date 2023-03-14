import os
from os import walk

import pandas as pd
import numpy as np
import sys
import re
import yaml
import argparse
import pickle

import numpy as np

from hla_typing_benchmark.parse_data import *

import http.client
from urllib.parse import urlparse


def checkUrl(url):
    p = urlparse(url)
    conn = http.client.HTTPConnection(p.netloc)
    conn.request('HEAD', p.path)
    resp = conn.getresponse()
    return resp.status < 400


def load_gs_data(gs_data_path = 'results/00_1000G_reference/1000G_hla_diversity_2014.txt', outfile_path = 'results/01_1000G_reference/1000G_2014_cleaned.tsv'):
    #Load gold standard data (refers to output from Snakemake run)
    #If Snakemake has not been run, the dataset can be found at http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20140725_hla_genotypes/20140702_hla_diversity.txt
    MG_exome_df = pd.read_csv(gs_data_path, sep = " ", comment='#')

    #Remove quotes
    #Change name of Utah individuals from CEPH to CEU as seen in the 1000 genomes database:
    #Remove samples with NaN (non typed alleles)
    MG_exome_df = MG_exome_df.replace({'\"':''}, regex=True).replace('CEPH','CEU').dropna()     

    #Only include Sample ID's that has been typed in the pipeline (read from config)
    with open("/work/nthu/publications/springers-hla-benchmark/snakemake/config.yaml", 'r') as stream:
        config = yaml.safe_load(stream)
    MG_exome_df = MG_exome_df[MG_exome_df['id'].isin(list(config['sample_urls']))].drop('sbgroup', axis = 1)

    #Replace 0000 with empty string:
    MG_exome_df = MG_exome_df.replace('0000', '')

    #Merge rows, where a person has been typed twice.
    #Check for non-identical rows
    print(f"Samples in gold standard dataset: {len(set(list(MG_exome_df['id'])))}")

    non_unique = list({x for x in list(MG_exome_df['id']) if list(MG_exome_df['id']).count(x) > 1})
    non_unique_df = MG_exome_df[MG_exome_df['id'].isin(non_unique)]

    clean_duplicates_df = pd.DataFrame()
    for column in ['A','A.1','B','B.1','C','C.1','DRB1','DRB1.1','DQB1','DQB1.1']:
        clean_duplicates_df[column] = non_unique_df.groupby(['id'])[column].apply('/'.join)
        
        #Remove potentaial starting '/'
        for identity in clean_duplicates_df.index:
            entry = clean_duplicates_df.loc[identity, column]
            
            if entry.startswith('/'):
                clean_duplicates_df.at[identity, column] = entry[1:]
            if entry.endswith('/'):
                clean_duplicates_df.at[identity, column] = entry[:-1]
            
    clean_duplicates_df = clean_duplicates_df.reset_index(drop=True)

    #Remove the duplicate rows from the full dataframe"
    MG_exome_df = MG_exome_df[~MG_exome_df['id'].isin(non_unique)]

    #Add back the clean duplicate rows:
    MG_exome_df = MG_exome_df.append(clean_duplicates_df, sort=False)

    #Reset index
    MG_exome_df = MG_exome_df.reset_index(drop = True)

    #Check that only uniwue entries exist now.
    assert len(set(list(MG_exome_df['id']))) == len(list(MG_exome_df['id']))

    #Set id as index
    MG_exome_df = MG_exome_df.set_index('id')


    #Look  for samples without gold standard data from 2014 Gorroud
    #Remember entries, which are not typed
    non_typed_samples = list()

    #Furthermore - in cases, where 2 field resolution is not available, use 2018 data.
    for identity in list(MG_exome_df.index):
        
        for col in MG_exome_df.columns:
            #get list of predictions
            predictions_list = MG_exome_df.loc[identity,col].split('/')
            
            #Put gene name in front of all entries in the list
            for i in range(len(predictions_list)):
                predictions_list[i] = col.split('.')[0] + "*" + predictions_list[i]
            
            #Check, that all entries in the list have at least four field resolution:
            for pred in predictions_list:
                
                #Remove all non valid entries
                if convert_to_two_field(pred) == None:
                    predictions_list.remove(pred)
            
            #Note, if no proper typing was made in 2014
            if len(predictions_list) == 0:
                non_typed_samples.append(identity)
                predictions_list = ['not_typed_in_2014']
                            
    
            MG_exome_df.at[identity,col] = predictions_list
    
    if len(non_typed_samples) > 0:
        print(f"Non typed samples in 2014: {non_typed_samples}")

    assert len(non_typed_samples) == 0

    #Merge haplotypes
    MG_exome_df['A_merged']= MG_exome_df[['A', 'A.1']].apply(lambda x: list(x), axis=1)
    MG_exome_df['B_merged']= MG_exome_df[['B', 'B.1']].apply(lambda x: list(x), axis=1)
    MG_exome_df['C_merged']= MG_exome_df[['C', 'C.1']].apply(lambda x: list(x), axis=1)
    MG_exome_df['DRB1_merged']= MG_exome_df[['DRB1', 'DRB1.1']].apply(lambda x: list(x), axis=1)
    MG_exome_df['DQB1_merged']= MG_exome_df[['DQB1', 'DQB1.1']].apply(lambda x: list(x), axis=1)

    MG_exome_merged_df = MG_exome_df.drop(columns=['A', 'A.1', 'B', 'B.1', 'C', 'C.1', 'DRB1', 'DRB1.1','DQB1', 'DQB1.1'])

    MG_exome_merged_df = MG_exome_merged_df.rename(columns={"A_merged": "A", "B_merged": "B", "C_merged": "C", "DRB1_merged": "DRB1", "DQB1_merged": "DQB1"})

    #Convert all predictions to 2-field resolution 
    gs_two_field_df = MG_exome_merged_df.copy()

    #Loop over all entries and update the three gold standard dataframes, so they fit the idividual resolution
    for identity in list(MG_exome_merged_df.index):
        for gene in list(MG_exome_merged_df.columns):
            
            old_pred_list = MG_exome_merged_df.loc[identity,gene]
            gene_pred_two_field = []
        
            #Loop over the two alleles
            for allele_list in old_pred_list:             
                allele_pred_two_field = []

                #Convert each prediction to it's respective correct format:
                for allele in allele_list:
                    allele_pred_two_field.append(convert_to_two_field(allele))

                #Merge the two alleles for each gene in a list
                gene_pred_two_field.append(list(set(allele_pred_two_field)))
            
            #Update dataframes wit0h the new predictions
            gs_two_field_df.at[identity,gene] = gene_pred_two_field

    gs_two_field_df.to_pickle(outfile_path)



def main():
    parser = get_argparser()
    args = parser.parse_args()
    load_gs_data(gs_data_path=args.input, outfile_path=args.output)


def get_argparser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--input',
                         help='Path to 1000G_hla_diversity_2014.txt. \nShould be downloaded first (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20140725_hla_genotypes/20140702_hla_diversity.txt)',
                         default='results/00_1000G_reference/1000G_hla_diversity_2014.txt',
                         required=False)

    parser.add_argument('--output',
                        help='Path to write formatted gold standard dataset.',
                        default='results/01_1000G_reference/1000G_2014_cleaned.pkl',
                        required=False)
    return parser


if __name__ == '__main__':
    main()