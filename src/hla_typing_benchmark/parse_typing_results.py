import os
import pandas as pd
import numpy as np
import sys
import re
import yaml
from yaml.loader import SafeLoader
import argparse
import json
import pickle
import matplotlib.pyplot as plt
from collections import Counter

from hla_typing_benchmark.parse_data import *


#Load resolution conversion dicts globally to avoid loading them every time convert_allele is called
p_group_dict = make_p_group_dict()
e_group_dict = make_e_group_dict()


def flatten(a):
    return [item for sublist in a for item in sublist if item != []]


def convert_allele(allele, resolution = 'two_field'):

    """
    input:
    allele (str): An allele in the format (A|B|C|DRB1|DQB1)\*\d{2}:\d{2,3}:?\d{0,3}G?:?\d{0,3} to be converted.
    resolution:   the resolution, with which the allele is converted to

    """
    if resolution == 'one_field':
        converted_allele = convert_to_one_field(allele)
    
    elif resolution == 'two_field':
        converted_allele = convert_to_two_field(allele)
        
    elif resolution == 'p_group':
        converted_allele = convert_to_p_group(allele, p_group_dict=p_group_dict)
    
    elif resolution == 'e_group':
        converted_allele = convert_to_e_group(allele, e_group_dict=e_group_dict, p_group_dict=p_group_dict)
            
    else:
        print('A conversion mistake happend. Please specify a correct conversion type.')
        converted_allele = None
        
    return converted_allele


def find_tool_results(tool_resultpath, result_extension = '_result.tsv'):
    tool_files = []
    for dirpath, subdirs, files in os.walk(tool_resultpath):
        for x in files:
            #Don't include the performance logs
            if x.endswith(result_extension):
                tool_files.append(os.path.join(dirpath, x))

    return tool_files


def load_kourami_results(kourami_result_filepath):

    kourami_files = find_tool_results(kourami_result_filepath)

    #Initalize result dict for single guess and for multiple typing
    kourami_results = {}

    for filename in kourami_files:
        subject_id = filename.split('/')[-1].replace('_result.tsv', '')

        #If file is not empty:
        if os.stat(filename).st_size != 0:
            temp_result_dict = {}    

            with open(filename, 'r') as infile:
                for line in infile:
                    #Find the first match / allele prediction
                    allele_searcher = re.search(r'(A|B|C|DRB1|DQB1)\*\d{2}:\d{2,3}:?\d{0,3}G?:?\d{0,3}', line)
                                            
                    if allele_searcher is not None:                        
                        found_allele = allele_searcher.group(0)
                        
                        gene = re.search(r'(A|B|C|DRB1|DQB1)', found_allele).group(0)
                        
                        #Convert alleles to right resolution:
                        pred_converted = convert_allele(found_allele)

                        if pred_converted != convert_to_two_field(found_allele):
                            print(filename)

                        #Add to list of predictions for this sample
                        if gene not in temp_result_dict.keys():
                            temp_result_dict[gene] = [[pred_converted]]
                        else:
                            temp_result_dict[gene] += [[pred_converted]]

            #Add sample prediction to dict
            kourami_results[subject_id] = temp_result_dict
            
        #If file is empty, add an empty dict.
        else:
            kourami_results[subject_id] = {}
            #print(f'Kourami: No prediction made for {filename} ')                           
                
    return kourami_results


def load_hla_la_results(hla_la_result_filepath):

    hla_la_files = find_tool_results(hla_la_result_filepath, 'R1_bestguess_G.txt')
    hla_la_results = {}

    for filename in hla_la_files:
        subject_id = filename.split('/')[2]
        temp_results_object = pd.read_csv(filename, sep = '\t')['Allele']
        temp_results = [i for i in temp_results_object if i.startswith(('A', 'B', 'C', 'DRB1', 'DQB1'))]
        
        #Make dict of dicts for results:
        temp_result_dict = {}
        
        for pred in temp_results:
            gene = re.search(r'(A|B|C|DRB1|DQB1)', pred).group(0)
            
            pred_converted = convert_allele(pred)
            #pred_converted = pred          
            
            #Add to list of predictions for this sample
            if gene not in temp_result_dict.keys():
                temp_result_dict[gene] = [[pred_converted]]
            else:
                temp_result_dict[gene] += [[pred_converted]]

        hla_la_results[subject_id] = temp_result_dict

    return hla_la_results


def load_optitype_results(optitype_result_filepath):

    optitype_files = find_tool_results(optitype_result_filepath)        
    optitype_results = {}

    for filename in optitype_files:
        subject_id = filename.split('/')[-1].replace('_result.tsv', '')
        
        #Check for right file, and that the file is not empty
        if os.stat(filename).st_size != 0:
            temp_results_raw = list(pd.read_csv(filename, sep = '\t').iloc[0])[1:7]

            temp_results = [i for i in temp_results_raw if isinstance(i,str)]
    
            #Make dict of dicts for results:
            temp_result_dict = {}
            
            for pred in temp_results:
                gene = re.search(r'(A|B|C|DRB1|DQB1)', pred).group(0)
                
                pred_converted = convert_allele(pred)
                #pred_converted = pred
                
                #Add to list of predictions for this sample
                if gene not in temp_result_dict.keys():
                    temp_result_dict[gene] = [[pred_converted]]
                else:
                    temp_result_dict[gene] += [[pred_converted]]
        
            optitype_results[subject_id] = temp_result_dict
    
        #Add empty entry for empty file
        else:
            optitype_results[subject_id] = {}
            
    return optitype_results


def load_hisat_genotype_results(hisatgenotype_result_filepath):

    hisatgenotype_files = find_tool_results(hisatgenotype_result_filepath, 'results.txt')
   
    #Save two predictions. One, with one guess per allele and one with the full prediction
    hisatgenotype_results = {}

    for filename in hisatgenotype_files:
        subject_id = filename.split('/')[2]
        hisatgenotype_resultlist = list()
            
        with open(filename) as infile:
            for line in infile:
                result = re.match(r'^\t+(1|2)\sranked (A|B|C|DRB1|DQB1)',line)
            
                if result is not None:
                    hisatgenotype_resultlist.append(line.split()[2])              
                    
            #Duplicate prediction for an allele in case of homologous case, so that each gene has two predictions.
            #In a homologous case, both result dicts only have one prediction and both needs an update.
            for allele in ['A', 'B', 'C', 'DRB1', 'DQB1']:
                allele_list = [pred for pred in hisatgenotype_resultlist if pred.startswith(allele)]

                if len(allele_list) == 1:
                    hisatgenotype_resultlist.append(allele_list[0])
                    hisatgenotype_resultlist.sort()

            temp_results = hisatgenotype_resultlist
            
        #Make dict of dicts for results:
        temp_result_dict = {}
        
        for pred in temp_results:
            gene = re.search(r'(A|B|C|DRB1|DQB1)', pred).group(0)
            
            pred_converted = convert_allele(pred)
            #pred_converted = pred
            
            #Add to list of predictions for this sample
            if gene not in temp_result_dict.keys():
                temp_result_dict[gene] = [[pred_converted]]
            else:
                temp_result_dict[gene] += [[pred_converted]]


        hisatgenotype_results[subject_id] = temp_result_dict
                    
    return hisatgenotype_results


def validate_call(correct_alleles, predicted_alleles, resolution):

    #Start by converting the alleles to the correct resolution
    correct_call_1 = {convert_allele(allele, resolution=resolution) for allele in correct_alleles[0]}
    correct_call_2 = {convert_allele(allele, resolution=resolution) for allele in correct_alleles[1]}
    correct_call = [list(correct_call_1), list(correct_call_2)]
    
    if predicted_alleles == '':
        num_correct_hits = 0
        pred_call = predicted_alleles
        return num_correct_hits, correct_call, pred_call

    pred_1 = {convert_allele(predicted_alleles[0][0], resolution=resolution)}
    pred_2 = {convert_allele(predicted_alleles[1][0], resolution=resolution)}
    pred_call = [list(pred_1), list(pred_2)]

    try:
        correct_0_pred_0 = list(correct_call_1.intersection(pred_1))
        correct_1_pred_1  = list(correct_call_2.intersection(pred_2))

        option_1 = [correct_0_pred_0, correct_1_pred_1]

        correct_0_pred_1  = list(correct_call_1.intersection(pred_2))
        correct_1_pred_0  = list(correct_call_2.intersection(pred_1))

        option_2 = [correct_0_pred_1, correct_1_pred_0]

        hits_1 = len([i for i in option_1 if i != []])
        hits_2 = len([i for i in option_2 if i != []])

        num_correct_hits = max(hits_1,hits_2)
        
    except KeyError as error:
        num_correct_hits = 0

    return num_correct_hits, correct_call, pred_call


#Get the number of predictions from a tool for a specific allele for a specific subject (2 or 0)
def get_count(tool_prediction, locus, subject):
    try:
        pred = tool_prediction[subject][locus]
        return len(pred)
    except KeyError as error:
        return 0


def validate_typing(typing_results_dict, gold_standard_df, subject_id_list, resolution):
    results_dict = {}
    full_typing_results_dict = {}

    #Loop over subjects, loci and tools
    for tool in typing_results_dict:
        results_dict[tool] = {}
        full_typing_results_dict[tool] = {}
        
        for locus in gold_standard_df.columns:
            results_dict[tool][locus] = {}
            results_dict[tool][locus]['count'] = 0
            results_dict[tool][locus]['score'] = 0

            full_typing_results_dict[tool][locus] = {}

            #Find counts of calls and correct calls
            for subject in subject_id_list:
                predicted_alleles = ''
                #If typing exists:
                if subject in typing_results_dict[tool]:
                    if locus in typing_results_dict[tool][subject]:
                        #Count number of total calls:
                        results_dict[tool][locus]['count'] += get_count(typing_results_dict[tool], locus, subject)
                    
                        #check whether it is valid.
                        predicted_alleles = typing_results_dict[tool][subject][locus]
                
                correct_alleles_list = gold_standard_df.loc[subject, locus]

                num_correct_hits, correct_call, pred_call = validate_call(correct_alleles_list, predicted_alleles, resolution)

                full_typing_results_dict[tool][locus][subject] = {}
                full_typing_results_dict[tool][locus][subject]['reference'] = correct_call
                full_typing_results_dict[tool][locus][subject]['prediction'] = pred_call
                full_typing_results_dict[tool][locus][subject]['miscalls'] = 2-num_correct_hits
                
                results_dict[tool][locus]['score'] += num_correct_hits


        #Add entries for all class I, class II and all alleles
        results_dict[tool]['HLA-I'] = {}
        results_dict[tool]['HLA-II'] = {}
        results_dict[tool]['Total'] = {}

        results_dict[tool]['HLA-I']['count'] = sum([results_dict[tool][locus]['count'] for locus in ['A', 'B', 'C']]) 
        results_dict[tool]['HLA-II']['count'] = sum([results_dict[tool][locus]['count'] for locus in ['DRB1', 'DQB1']]) 
        results_dict[tool]['Total']['count'] = sum([results_dict[tool][locus]['count'] for locus in ['HLA-I', 'HLA-II']]) 

        results_dict[tool]['HLA-I']['score'] = sum([results_dict[tool][locus]['score'] for locus in ['A', 'B', 'C']]) 
        results_dict[tool]['HLA-II']['score'] = sum([results_dict[tool][locus]['score'] for locus in ['DRB1', 'DQB1']]) 
        results_dict[tool]['Total']['score'] = sum([results_dict[tool][locus]['score'] for locus in ['HLA-I', 'HLA-II']]) 


        #Find total counts and calculate call rate and typing accuracy
        for locus in results_dict[tool]:
            if locus in ['A', 'B', 'C', 'DRB1', 'DQB1']:
                multiplier = 1
            elif locus == 'HLA-I':
                multiplier = 3
            elif locus == 'HLA-II':
                multiplier = 2
            elif locus == 'Total':
                multiplier = 5

            total_correct_calls = multiplier * len(gold_standard_df[gold_standard_df.index.isin(subject_id_list)]) * 2
            results_dict[tool][locus]['call_rate'] = results_dict[tool][locus]['count'] * 100  / total_correct_calls
            results_dict[tool][locus]['typing_accuracy'] = results_dict[tool][locus]['score'] * 100 / total_correct_calls


    return results_dict, full_typing_results_dict


def load_all_results(gs_data = 'results/01_1000G_reference/1000G_2014_cleaned.pkl',
                        kourami_path = None, 
                        hla_la_path = None,
                        optitype_path = None, 
                        hisat_genotype_path = None):

    gs_two_field_df = pd.read_pickle(gs_data)

    typing_results_dict = {}


    #Load results of the tools:
    if kourami_path != None:
        typing_results_dict['Kourami'] = load_kourami_results(kourami_path)
    
    if hla_la_path != None:
        typing_results_dict['HLA-LA'] = load_hla_la_results(hla_la_path)
    
    if optitype_path != None:
        typing_results_dict['Optitype'] = load_optitype_results(optitype_path)
    
    if hisat_genotype_path != None:
        typing_results_dict['Hisatgenotype'] = load_hisat_genotype_results(hisat_genotype_path)

    rename_resolutions = {
        'one_field' : '1-field',
        'e_group' : 'pseudosequence',
        'p_group' : 'P group',
        'two_field' : '2-field',
    }

    #Load list of gold standard samples throug the config
    with open('snakemake/config.yaml', 'r') as f:
        config = yaml.load(f, Loader=SafeLoader)
    
    gold_standard_id_list = list(config['sample_urls'].keys())

    #Loop over typing resolutions 
    four_resolutions_results = {}
    full_typing_results_dict = {}

    for resolution in ['one_field', 'e_group', 'p_group', 'two_field']:
        four_resolutions_results[rename_resolutions[resolution]], full_typing_results_dict[rename_resolutions[resolution]] = validate_typing(typing_results_dict=typing_results_dict, gold_standard_df=gs_two_field_df, subject_id_list=gold_standard_id_list, resolution=resolution)

    return four_resolutions_results
 
 


def main():
    parser = get_argparser()
    args = parser.parse_args()

    full_results = load_all_results(gs_data = args.gs_data,
                                    kourami_path=args.kourami,
                                    hla_la_path=args.hla_la,
                                    optitype_path=args.optitype,
                                    hisat_genotype_path=args.hisat_genotype)


    if not os.path.exists(os.path.dirname(args.output)):
        os.makedirs(os.path.dirname(args.output))
        
    with open(args.output, 'w', encoding='utf-8') as f:
        json.dump(full_results, f, ensure_ascii=False, indent=4)



def get_argparser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--gs-data',
                        help='Gold standard data (.pkl Pandas dataframe file)',
                        default='results/01_1000G_reference/1000G_2014_cleaned.pkl')


    parser.add_argument('--kourami',
                         help='Folder with HLA typing results from Kourami')

    parser.add_argument('--optitype',
                         help='Folder with HLA typing results from Optitype')

    parser.add_argument('--hla-la',
                         help='Folder with HLA typing results from HLA*LA')

    parser.add_argument('--hisat-genotype',
                         help='Folder with HLA typing results from HISAT-genotype')

    parser.add_argument('--output',
                        help='Path to save .json file with complete, collected typing results including call rate, typing accuracy etc.',
                        default='results/05_collected_typing_results/typing_results.json')
    return parser


if __name__ == '__main__':
    main()