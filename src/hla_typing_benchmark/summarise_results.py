import numpy as np
import sys
import re
import json
import math
import pandas as pd
idx = pd.IndexSlice

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
import seaborn as sns
import argparse

from collections import Counter


#Set seaborn color palette:
palette_list = list(sns.color_palette("colorblind"))


def generate_results_as_table(results_dict):

    #Reformat json structure into pandas Dataframe
    reformatted_dict = dict()
    for resolution, tool_dict in results_dict.items():
        for tool, locus_dict in tool_dict.items():
            for locus, metric_dict in locus_dict.items():
                for metric, value in metric_dict.items():
                    if (tool, resolution) not in reformatted_dict:
                        reformatted_dict[(tool, resolution)] = dict()

                    #Drop count, score and save typing accuracy per resolution
                    if metric == 'typing_accuracy':
                        reformatted_dict[(tool, resolution)][locus] = value
                    elif (metric == 'call_rate'):
                        if (tool, metric) not in reformatted_dict:
                            reformatted_dict[(tool, metric)] = dict()
                        if locus not in reformatted_dict[(tool, metric)]:
                            reformatted_dict[(tool, metric)][locus] = value

        results_df = pd.DataFrame(reformatted_dict)


        tool_order = ['Kourami', 'HLA-LA', 'Optitype', 'Hisatgenotype', 'STC-seq']
        metric_order = ['call_rate', '1-field', 'pseudosequence', 'P group', '2-field']

        multi_tuples = []
        for tool in tool_order:
            for metric in metric_order:
                multi_tuples += [(tool, metric)]

        multi_cols = pd.MultiIndex.from_tuples(multi_tuples, names=['Tool', 'Metric'])

        results_df = pd.DataFrame(results_df, columns=multi_cols)

    return results_df


def extract_results(results_dict, allele_index, resolution, labels):
    res_list = [results_dict[resolution][l][allele_index]['typing_accuracy'] for l in labels]
    return res_list 


# Plots of the Performance of the tools

def make_class_plot_from_allele_list(results_dict):

    labels1 = ['Kourami', 'HLA-LA', 'Hisatgenotype', 'Optitype']
    labels2 = ['Kourami', 'HLA-LA', 'Hisatgenotype',]
    textlabels1 = ['Kourami', 'HLA*LA', 'HISAT-genotype', 'Optitype']
    textlabels2 = ['Kourami', 'HLA*LA', 'HISAT-genotype', '']

    fig, axs = plt.subplot_mosaic([['A'], ['B']], constrained_layout=True, figsize=(18,12))

    for label, ax in axs.items():

        if label == 'A':
            allele_index = 'HLA-I'
            title = 'HLA Class I alleles (HLA-A, -B and -C)'
        else:
            allele_index = 'HLA-II'
            title = 'HLA Class II alleles (HLA-DRB1 and -DQB1)'

        call_rate = [results_dict['1-field'][l][allele_index]['call_rate'] for l in labels1]
        accuracy_one_field =  extract_results(results_dict, allele_index, '1-field', labels1)
        accuracy_two_field = extract_results(results_dict, allele_index, '2-field', labels1)
        accuracy_p_group = extract_results(results_dict, allele_index, 'P group', labels1)
        accuracy_e_group =  extract_results(results_dict, allele_index, 'pseudosequence', labels1)
    
        x = np.arange(len(labels1)) # the label locations
        
    
        width = 0.18  # the width of the bars
        rects1 = ax.bar(x - 10.3*width/5, call_rate, width, label='call rate', color = '#808080')
        rects2 = ax.bar(x - 5*width/5, accuracy_one_field, width, label='1-field accuracy', color = palette_list[0])
        rects3 = ax.bar(x, accuracy_e_group, width, label='Pseudoseq accuracy', color = palette_list[1])
        rects4 = ax.bar(x + 5*width/5, accuracy_p_group, width, label='P group accuracy', color = palette_list[2])
        rects5 = ax.bar(x + 10*width/5, accuracy_two_field, width, label='2-field accuracy', color = palette_list[3])

        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax.set_ylabel('%', size = 22)
        ax.set_xticks(x)
        
        if allele_index in ('HLA-II', 'DRB1', 'DQB1'):
            ax.set_xticklabels(textlabels2, size = 18)
        else:
            ax.set_xticklabels(textlabels1, size = 18)

        ax.legend()


        def autolabel(rects, y_height):
            """Attach a text label above each bar in *rects*, displaying its height."""
            for rect in rects:
                if (allele_index in ('HLA-II', 'DRB1', 'DQB1')) and (rect.get_x() > 2.5):
                    continue
                else:
                    height = rect.get_height()
                    ax.annotate(f"{format(height, '#.3g')}",
                                xy=(rect.get_x() + rect.get_width() / 2, y_height),
                                xytext=(0, 3),  # 3 points vertical offset
                                textcoords="offset points",
                                ha='center', va='bottom', 
                                color="black", size = 16, weight =400)

        autolabel(rects1, 50)
        autolabel(rects2, 40)
        autolabel(rects3, 30)
        autolabel(rects4, 20)
        autolabel(rects5, 10)        
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)       
        
        ax.grid(which='minor', alpha=0.4)
        # Major ticks every 20, minor ticks every 5
        minor_ticks = np.arange(0, 101, 5)

        ax.set_yticks(minor_ticks, minor=True)
        ax.tick_params(axis="y", labelsize=16)


        # Or if you want different settings for the grids:
        ax.grid(axis = 'y')
        ax.set_axisbelow(True)

        if allele_index in ('HLA-II', 'DRB1', 'DQB1'):
            ax.legend(loc='right', prop={'size': 18})
        else:
            ax.get_legend().remove()
        plt.tight_layout()

        #plt.title('HLA typing performance for HLA-A, -B, -C, -DRB1 and DQB1', size = 24)
        #plt.title('HLA typing performance for HLA class II genes (with ensemble)', size = 24)
        ax.set_title(f'HLA typing performance for {title}', size = 22)
        trans = mtransforms.ScaledTranslation(-20/72, 7/72, fig.dpi_scale_trans)
        ax.text(-0.01, 1.0, label, transform=ax.transAxes + trans,
                fontsize=22, va='bottom', fontfamily='serif', weight='bold')
    
    plt.tight_layout()
    
    return plt





def make_loci_plot_from_allele_list(results_dict):

    labels1 = ['Kourami', 'HLA-LA', 'Hisatgenotype', 'Optitype']
    labels2 = ['Kourami', 'HLA-LA', 'Hisatgenotype', 'Optitype']
    textlabels1 = ['Kourami', 'HLA*LA', 'HISAT-genotype', 'Optitype']
    textlabels2 = ['Kourami', 'HLA*LA', 'HISAT-genotype', '']

    fig, axs = plt.subplot_mosaic([['A'], ['B'], ['C'], ['D'], ['E']], constrained_layout=True, figsize=(18,24))

    for label, ax in axs.items():

        if label == 'A':
            allele_index = 'A'
            title = 'HLA-A'
        elif label == 'B':
            allele_index = 'B'
            title = 'HLA-B'
        elif label == 'C':
            allele_index = 'C'
            title = 'HLA-C'
        elif label == 'D':
            allele_index = 'DRB1'
            title = 'HLA-DRB1'
        elif label == 'E':
            allele_index = 'DQB1'
            title = 'HLA-DQB1'

        if allele_index in ('HLA-II', 'DRB1', 'DQB1'):
            call_rate = [results_dict['1-field'][l][allele_index]['call_rate'] for l in labels2]
            accuracy_one_field =  extract_results(results_dict, allele_index, '1-field', labels2)
            accuracy_two_field = extract_results(results_dict, allele_index, '2-field', labels2)
            accuracy_p_group = extract_results(results_dict, allele_index, 'P group', labels2)
            accuracy_e_group =  extract_results(results_dict, allele_index, 'pseudosequence', labels2)
        
            x = np.arange(len(labels2))
            
        else:   
            call_rate = [results_dict['1-field'][l][allele_index]['call_rate'] for l in labels1]
            accuracy_one_field =  extract_results(results_dict, allele_index, '1-field', labels1)
            accuracy_two_field = extract_results(results_dict, allele_index, '2-field', labels1)
            accuracy_p_group = extract_results(results_dict, allele_index, 'P group', labels1)
            accuracy_e_group =  extract_results(results_dict, allele_index, 'pseudosequence', labels1)

            x = np.arange(len(labels1))  # the label locations
    

        width = 0.18  # the width of the bars
        rects1 = ax.bar(x - 10.4*width/5, call_rate, width, label='call rate', color = '#808080')
        rects2 = ax.bar(x - 5*width/5, accuracy_one_field, width, label='1-field accuracy', color = palette_list[0])
        rects3 = ax.bar(x, accuracy_e_group, width, label='Pseudoseq accuracy', color = palette_list[1])
        rects4 = ax.bar(x + 5*width/5, accuracy_p_group, width, label='P group accuracy', color = palette_list[2])
        rects5 = ax.bar(x + 10*width/5, accuracy_two_field, width, label='2-field accuracy', color = palette_list[3])

        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax.set_ylabel('%', size = 20)
        ax.set_xticks(x)
        
        if allele_index in ('HLA-II', 'DRB1', 'DQB1'):
            ax.set_xticklabels(textlabels2, size = 18)
        else:
            ax.set_xticklabels(textlabels1, size = 18)

        ax.legend()


        def autolabel(rects, y_height):
            """Attach a text label above each bar in *rects*, displaying its height."""
            for rect in rects:
                if (allele_index in ('HLA-II', 'DRB1', 'DQB1')) and (rect.get_x() > 2.5):
                    continue
                else:
                    height = rect.get_height()
                    ax.annotate(f"{format(height, '#.3g')}",
                                xy=(rect.get_x() + rect.get_width() / 2, y_height),
                                xytext=(0, 3),  # 3 points vertical offset
                                textcoords="offset points",
                                ha='center', va='bottom', 
                                color="black", size = 16, weight =400)

        autolabel(rects1, 50)
        autolabel(rects2, 40)
        autolabel(rects3, 30)
        autolabel(rects4, 20)
        autolabel(rects5, 10)    
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)       
        
        ax.grid(which='minor', alpha=0.4)

        minor_ticks = np.arange(0, 101, 5)

        ax.set_yticks(minor_ticks, minor=True)
        ax.tick_params(axis="y", labelsize=16)


        # Or if you want different settings for the grids:
        ax.grid(axis = 'y')
        ax.set_axisbelow(True)

        if allele_index == 'DQB1':
            ax.legend(loc='lower right', prop={'size': 16})
        else:
            ax.get_legend().remove()
        plt.tight_layout()

        ax.set_title(f'HLA typing performance for {title}', size = 22)
        trans = mtransforms.ScaledTranslation(-20/72, 7/72, fig.dpi_scale_trans)
        ax.text(-0.01, 1.0, label, transform=ax.transAxes + trans,
                fontsize=22, va='bottom', fontfamily='serif', weight='bold')
    

    plt.tight_layout()
    
    return plt





def main():
    parser = get_argparser()
    args = parser.parse_args()

    #Load and reformat results
    with open(args.typing_results, 'r') as infile:
        results_dict = json.load(infile)

    results_df = generate_results_as_table(results_dict)

    tool_order = ['Kourami', 'HLA-LA', 'Optitype', 'Hisatgenotype']
    metric_order = ['call_rate', '1-field', 'pseudosequence', 'P group', '2-field']

    multi_tuples = []
    for tool in tool_order:
        for metric in metric_order:
            multi_tuples += [(tool, metric)]

    multi_cols = pd.MultiIndex.from_tuples(multi_tuples, names=['Tool', 'Metric'])

    results_df = pd.DataFrame(results_df, columns=multi_cols)

    results_df.to_csv(args.results_table, sep='\t')


    #Generate plots

    fig = make_class_plot_from_allele_list(results_dict)
    fig.savefig(args.class_plot, dpi=600)

    fig = make_loci_plot_from_allele_list(results_dict)
    fig.savefig(args.loci_plot, dpi=600)



def get_argparser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--typing-results',
                         help='Path to typing_results.json (output from parse_typing_results)',
                         default='results/05_collected_typing_results/typing_results.json')

    parser.add_argument('--results-table',
                        help='Path to write a table with a summary of results',
                        default='results/01_1000G_reference/1000G_2014_cleaned.tsv')

    parser.add_argument('--class-plot',
                        help='Generate plot with an overview of HLA Class I and HLA class II performance',
                        default='results/05_collected_typing_results/HLA_class_performance.jpg')

    parser.add_argument('--loci-plot',
                        help='Generate plot with an overview of HLA-A, -B, -C and -DRB1 performance',
                        default='results/05_collected_typing_results/HLA_loci_performance.jpg')

    return parser


if __name__ == '__main__':
    main()