import http.client
from urllib.parse import urlparse


def checkUrl(url):
    p = urlparse(url)
    conn = http.client.HTTPConnection(p.netloc)
    conn.request('HEAD', p.path)
    resp = conn.getresponse()
    return resp.status < 400


#Check all exome links to create a list of id's with exome data
gold_standard_id_list = []
gold_standard_url_list = []

#Load gold standard data (refers to output from Snakemake run)
#If Snakemake has not been run, the dataset can be found at http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20140725_hla_genotypes/20140702_hla_diversity.txt
gs_2014_df = pd.read_csv('../../results/00_1000G_reference/1000G_hla_diversity_2014.txt', sep = " ", comment='#')

for sample_id in gs_2014_df['id']:
    sbgroup = gs_2014_df[gs_2014_df['id'] == sample_id]['sbgroup'].iloc[0]         
    wget_url = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/{}/{}/exome_alignment/{}.alt_bwamem_GRCh38DH.20150826.{}.exome.cram".format(sbgroup, identity, identity, sbgroup, )
    
    if checkUrl(wget_url):
        gold_standard_id_list.append(sample_id)
        gold_standard_url_list.append(wget_url)


with open(configfile_1, 'w') as outfile:
    outfile.write('sample_urls:\n')
    for entry in full_sample_urllist:
        entry_name = entry.split('/')[9]
        outfile.write("  " + entry_name + ": " +  "\"" + entry + "\"" + '\n')
    
    outfile.write('reference_genome:\n\tGRCh38: "ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"')


