#!/usr/bin/env python3

import jinja2
import os
import argparse
import fnmatch
import pandas as pd
from urllib.parse import urljoin
import get_cosmic

############################ Variants paired tumor:normal samples from annotations.paired.tsv file ###########################

### Find the annotations paired file and its path, provided the run id and its path ###
def find(pattern, path):
    exclude_directories = set(['work'])
    file_path = []
    for root,dirs,files in os.walk(path):
        dirs[:] = [d for d in dirs if d not in exclude_directories]
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                file_path.append(os.path.join(root, name))
    return file_path

def sample_name(annotation_file,column):
    TS_number = annotation_file[column].str.split('_', expand=True,n=5)
    return TS_number[5]
 
def file_processing(annotations_paired_file):
    annotations_paired_file = annotations_paired_file.drop(columns=['NORMAL.DP','TUMOR.DP','TUMOR.AF','TUMOR.AD.REF','TUMOR.AD.ALT','TUMOR.AD.TOTAL','NORMAL.AD.REF','NORMAL.AD.ALT','NORMAL.AD.TOTAL','Run','Time','Session','Workflow','Location','System','GitBranch','GitTag'],axis = 1)
    annotations_paired_file['NORMAL.AF'].replace({'.': None}, inplace=True)
    annotations_paired_file = annotations_paired_file.astype({'NORMAL.AF':float, 'POS':str})
    annotations_paired_file = annotations_paired_file.round({'NORMAL.AF':3, 'AF':3})
    annotations_filtered_NAF = annotations_paired_file[(annotations_paired_file['NORMAL.AF'].isna()) | (annotations_paired_file['NORMAL.AF'] < 0.02)]
    annotations_processed = annotations_filtered_NAF.copy()
    annotations_processed['Variant'] = annotations_processed[['CHROM','POS','REF','ALT']].apply(lambda x: ':'.join(x),axis=1)
    annotations_processed['Bam_Variant'] = annotations_processed[['CHROM','POS']].apply(lambda x: '_'.join(x), axis=1)
    annotations_processed['MuTect2'] = annotations_processed['VariantCaller'].apply(lambda x: 'Yes' if x == "MuTect2" else 'No')
    annotations_processed['Strelka'] = annotations_processed['VariantCaller'].apply(lambda x: 'Yes' if x == "Strelka" else 'No')
    annotations_processed['LoFreqSomatic'] = annotations_processed['VariantCaller'].apply(lambda x: 'Yes' if x == "LoFreqSomatic" else 'No')
    return annotations_processed
    
def somatic_annotation_filtered(annotations_processed):
    final_df = annotations_processed
    final_dfs_nodup = final_df.drop_duplicates()
    sample_list = annotations_processed.Tumor.unique()

    variants = pd.DataFrame(columns=['Tumor','Normal','Gene.refGene','MuTect2','Strelka','LoFreqSomatic','AAChange.refGene','Variant','DP','AF','Func.refGene',
                          'ExonicFunc.refGene','NORMAL.AF'])
    cols = variants.columns.tolist()

    for sample in sample_list:
    ## Get all variants per sample ##
        test = final_dfs_nodup[(final_dfs_nodup['Tumor'] == sample)]
        final_df = test[cols]
        all_variants = pd.DataFrame(final_df)
        variants = pd.concat([variants,all_variants])
        variants.reset_index(drop=True, inplace=True)
    
    all_variants_final = variants.sort_values('Gene.refGene')
    all_variants_final_filter = all_variants_final[(all_variants_final['AF'] >= 0.05) & (all_variants_final['DP'] > 50)]
    all_variants_final_df = all_variants_final_filter[(all_variants_final_filter['ExonicFunc.refGene'] != 'synonymous SNV' ) & (all_variants_final_filter['Func.refGene'] == 'exonic') | (all_variants_final_filter['Func.refGene'] == 'upstream')| (all_variants_final_filter['Func.refGene'] == 'splicing')]
    return all_variants_final_df

def germline_lofreq_annotation(annotations_paired_file,lofreq_annotation_file,germline_genes_file):
    ### Consider only normal samples for germline filter ###
    normal_samples = annotations_paired_file.Normal.unique()
    normal_only_annotations_lofreq = pd.DataFrame(columns=['Sample','VariantCaller','Gene.refGene','CHROM','POS','REF','ALT','DP','AF','AAChange.refGene','CLINSIG','Func.refGene','ExonicFunc.refGene'])
    cols_lofreq = normal_only_annotations_lofreq.columns.tolist()

    for n in normal_samples:
        lofreq_annotation_file_normals = lofreq_annotation_file[lofreq_annotation_file['Sample'].str.contains(n)]
        normal_only_annotations_lofreq = pd.concat([normal_only_annotations_lofreq,lofreq_annotation_file_normals],sort=True)

    normal_only_annotations_lofreq['Sample'] = sample_name(normal_only_annotations_lofreq,'Sample')

    ### Check for only 59 genes, apply VAF > 30% and CLINSIG should contain pathogenic ###
    germline_genes = germline_genes_file['Germline_Gene'].values.tolist()
    germline_genes_lofreqannotation = normal_only_annotations_lofreq[normal_only_annotations_lofreq['Gene.refGene'].isin(germline_genes)]
    germline_genes_lofreqannotation_df = germline_genes_lofreqannotation.copy()
    germline_genes_lofreqannotation_filtered = germline_genes_lofreqannotation_df[(germline_genes_lofreqannotation_df['AF'] > 0.30)]
    germline_genes_lofreqannotation_filtered_df = germline_genes_lofreqannotation_filtered.copy()
    germline_genes_lofreqannotation_filtered_df['AF'] = germline_genes_lofreqannotation_filtered['AF'].round(decimals=3)
    germline_genes_lofreqannotation_filtered_clinsig = germline_genes_lofreqannotation_filtered_df[(germline_genes_lofreqannotation_filtered_df['CLINSIG'].str.contains('pathogenic')) | (germline_genes_lofreqannotation_filtered_df['CLINSIG'].str.contains('Likely pathogenic')) | (germline_genes_lofreqannotation_filtered_df['CLINSIG'].str.contains('Pathogenic')) | (germline_genes_lofreqannotation_filtered_df['CLINSIG'] == '.') & (germline_genes_lofreqannotation_filtered_df['Func.refGene'] == 'exonic') | (germline_genes_lofreqannotation_filtered_df['Func.refGene'] == 'upstream')| (germline_genes_lofreqannotation_filtered_df['Func.refGene'] == 'splicing')]

    germline_genes_lofreqannotation_filtered_clinsig_final = germline_genes_lofreqannotation_filtered_clinsig[(germline_genes_lofreqannotation_filtered_clinsig['ExonicFunc.refGene'] != 'synonymous SNV')]

    lofreq_germline_df = germline_genes_lofreqannotation_filtered_clinsig_final[cols_lofreq]
    return lofreq_germline_df
    
       
def main(runid,rundir,pactid):
    ## 1. process somatic variants ##
    somatic_annotation_file = find('annotations.paired.tsv', rundir)
    annotations_paired_file = pd.read_csv(somatic_annotation_file[0],sep="\t",header=0)
    annotations_processed = file_processing(annotations_paired_file)
    somatic_variants = somatic_annotation_filtered(annotations_processed)
    print("Somatic variants are generated!!")
    
     ## merge the ngs test case ID with the variants list ##
    maindir = rundir.rpartition('output')[0]
    demux_ss_find = find('demux-samplesheet.csv', maindir)
    demux_ss_file = pd.read_csv(demux_ss_find[0],skiprows=19)
    somatic_demux_merge = pd.merge(somatic_variants,demux_ss_file,left_on="Tumor",right_on="Sample_ID")
    somatic_demux_merge_excludedHAPMAP = somatic_demux_merge.loc[(somatic_demux_merge['Test_Number'] != '0')]
    somatic_demux_merge_final = somatic_demux_merge_excludedHAPMAP.copy()
    somatic_demux_merge_final['Tumor'] = sample_name(somatic_demux_merge_final,'Tumor')
    somatic_demux_merge_final['Normal'] = sample_name(somatic_demux_merge_final,'Normal')
    somatic_variants_final = somatic_demux_merge_final[['Test_Number','Tumor','Normal','Gene.refGene','MuTect2','Strelka','LoFreqSomatic','ExonicFunc.refGene','AAChange.refGene','Variant','DP','AF','Func.refGene','NORMAL.AF']]
    
    ## find per caller annotations to add cosmic info to report ##
    files_to_find = ['annotations.MuTect2.tsv', 'annotations.Strelka.tsv', 'annotations.LoFreqSomatic.tsv']
    files_per_caller = [find(file, rundir) for file in files_to_find]
    cosmic_results = [ get_cosmic.main(caller_file[0]) for caller_file in files_per_caller]
    cosmic_merged = pd.concat([cosmic_results[0], cosmic_results[1], cosmic_results[2]]).drop_duplicates()
    ## combine somatic variants df and cosmic variants df ##
    cosmic_somatic_merge = somatic_variants_final.merge(cosmic_merged,on = 'Variant', how = 'left').fillna(0)
    cosmic_somatic_merge['COSMIC_Count'] = cosmic_somatic_merge['COSMIC_Count'].astype(int)
    cosmic_somatic_merge['COSMIC_Count'] = cosmic_somatic_merge['COSMIC_Count'].astype(str).replace('0','None')
    cosmic_somatic_merge['CosmicID'] = cosmic_somatic_merge['CosmicID'].replace(0,'None')
    ## Get the required cols for final report ##
    somatic_cosmic_variants = cosmic_somatic_merge[['Test_Number','Tumor','Normal','Gene.refGene','Variant','DP','AF','MuTect2','Strelka','LoFreqSomatic','ExonicFunc.refGene','CosmicID','COSMIC_Count','AAChange.refGene','Func.refGene','NORMAL.AF']]
    somatic_file_out = '%s/%s/%s'% (rundir,"clinical","somatic_variants.csv")
    somatic_cosmic_variants.to_csv(somatic_file_out, index=False)

    somatic_variants_samples = somatic_variants_final.Test_Number.unique()

    ## 2. Process germline variants ##
    lofreq_annotation_file_find = find('annotations.LoFreq.tsv', rundir)
    lofreq_annotation_file = pd.read_csv(lofreq_annotation_file_find[0],sep="\t",header=0)
    germline_genes_file_find = find('germline_genes.csv', maindir)
    germline_genes_file = pd.read_csv(germline_genes_file_find[0],sep=",",header=0)
    germline_variants = germline_lofreq_annotation(annotations_paired_file,lofreq_annotation_file,germline_genes_file)
    print("Germline variants are generated!!")


#### Using Jinja2 to render the variants html template report ####

    template_path = os.path.join(maindir,"template/")
    templateLoader = jinja2.FileSystemLoader(searchpath=template_path)
    html_out_path = '%s/%s/'% (rundir,"clinical")
    templateEnv = jinja2.Environment(loader=templateLoader)
    TEMPLATE_FILE_SOMATIC = "somatic_variants.html"
    TEMPLATE_FILE_GERMLINE = "germline_variants.html"
    template_somatic = templateEnv.get_template(TEMPLATE_FILE_SOMATIC)
    template_germline = templateEnv.get_template(TEMPLATE_FILE_GERMLINE)

    html_out_somatic = template_somatic.render(somatic_variants=somatic_cosmic_variants,somatic_variants_samples=somatic_variants_samples,pactid=pactid,runid=runid)
    html_out_germline = template_germline.render(germline_variants=germline_variants,pactid=pactid,runid=runid)
    if not os.path.exists(html_out_path):
        os.mkdir(html_out_path)
    somatic_filename = '%s/%s-%s'% (html_out_path,pactid,"Somatic_Variants.html")
    html_file_somatic = open(somatic_filename, 'w')
    html_file_somatic.write(html_out_somatic)
    html_file_somatic.close()

    germline_filename = '%s/%s-%s'% (html_out_path,pactid,"Germline_Variants.html")
    html_file_germline = open(germline_filename, 'w')
    html_file_germline.write(html_out_germline)
    html_file_germline.close()
    os.system("chmod -R g+rwx " + html_out_path )
    print("Variants QC reports are rendered!!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Provide path to run directory's output dir")
    parser.add_argument('-rid','--runid', help='Your run id',required=True)
    parser.add_argument('-rdir','--rundir', help='Your run directory',required=True)
    parser.add_argument('-pactid','--pactid',help='PACT ID',required=True)
    args = parser.parse_args()
    main(args.runid,args.rundir,args.pactid)
