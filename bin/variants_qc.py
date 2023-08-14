#!/usr/bin/env python3

import re
import jinja2
import os,fnmatch,argparse,subprocess
import pandas as pd
from urllib.parse import urljoin

######################### Moving and processing IGV-snapshots dir #############################################

### Create dirs and move dirs to per caller type ###
def create_percaller_dirs(runid):
    dest_dir = "/gpfs/data/sequence/results/external/NYU/snuderllab/igv-snapshots/"
    full_run_path = '%s%s'% (dest_dir,runid)
    print(full_run_path)
    ## make run dir in minio browser ##
    if not os.path.exists(full_run_path):
        os.mkdir(full_run_path)
    each_caller_path = ["MuTect2","Strelka","LoFreqSomatic"] 
    caller_paths = []
    for caller in each_caller_path:
        caller_path = '%s/%s/'% (full_run_path,caller)
        if not os.path.exists(caller_path):
            os.mkdir(caller_path)
        caller_paths.append(caller_path)
    print("Done creating dirs!!")
    return caller_paths, full_run_path

    
def subprocess_cmd(command):
    p = subprocess.Popen(command,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    p.wait()

def move_dirs_percaller(rundir):
    igv_snapshotdir = '%s/%s/'% (rundir,"igv-snapshots")
    os.chdir(igv_snapshotdir)
    find_move_Dirs = "find . -maxdepth 1 -type d -name '*MuTect2*' -print0 | xargs -0 -I {} mv {} MuTect2/ | find . -maxdepth 1 -type d -name '*Strelka*' -print0 | xargs -0 -I {} mv {} Strelka/ | find . -maxdepth 1 -type d -name '*LoFreqSomatic*' -print0 | xargs -0 -I {} mv {} LoFreqSomatic/ "
    subprocess_cmd(find_move_Dirs)
    print("Moved per caller dirs!!")
    
def copy_png_toexternal(rundir,runid):
    print("QC Script started!!")
    igv_snapshotdir = '%s/%s/'% (rundir,"igv-snapshots")
    move_dirs_percaller(rundir)
    #caller_paths, full_run_path = create_percaller_dirs(runid)
    #os.chdir(igv_snapshotdir)
    #print("Copying png files to external results dir....")
    #find_png_copy_mutect2 = "find MuTect2/ -type f -name '*.png' -exec cp -R {} " + caller_paths[0] + " \;" 
    #find_png_copy_strelka = "find Strelka/ -type f -name '*.png' -exec cp -R {} " + caller_paths[1] + " \;"
    #find_png_copy_lofreqsomatic = "find LoFreqSomatic/ -type f -name '*.png' -exec cp -R {} " + caller_paths[2] + " \;"
    #subprocess_cmd(find_png_copy_mutect2)
    #subprocess_cmd(find_png_copy_strelka)
    #subprocess_cmd(find_png_copy_lofreqsomatic)
    #os.system("chmod -R g+rwx " + full_run_path )
    #print("Copied all png files to external dir!!")

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
    
def make_clickable(link,name):
    # target _blank to open new window
    return '<a target="_blank" href="{}">{}</a>'.format(link, name)

def igv_png(runID,igvsnapshot_path,callerType):
    png_caller = '%s_%s'% ('png',callerType)
    igvsnapshot_path_caller = '%s%s%s'% (igvsnapshot_path,"/igv-snapshots/",callerType)
    png_caller = []
    for root, dirs, files in os.walk(igvsnapshot_path_caller):
        for file in files:
            if file.endswith('.png'):
                png_caller.append(file)
    return png_caller

def add_igvsnapshots_to_variants(runID,igvsnapshot_path,callerType,annotations_processed,urlpath):
    png_caller = igv_png(runID,igvsnapshot_path,callerType)
    caller_variants_igv = pd.DataFrame(columns=['Tumor','Normal','Gene.refGene','AAChange.refGene','Variant','MuTect2','Strelka','LoFreqSomatic','DP','AF','Func.refGene','ExonicFunc.refGene','NORMAL.AF'])
    cols_final = caller_variants_igv.columns.tolist()
    for i in png_caller:
        groups = i.split('_')
        ## split the filename at 2nd occurence of '_' to match df['Bam_Variant'] column ##
        file_variant = '_'.join(groups[:2]), '_'.join(groups[2:])
        variants = annotations_processed[(annotations_processed[callerType] == 'Yes')]
        df = variants[variants['Bam_Variant'].str.contains(file_variant[0])]
        df_final = df.copy()
        urlpath_final = '%s/%s/%s/'% (urlpath,runID,callerType)
        url = urljoin(urlpath_final, i)
        df_final = df_final.drop(['Bam_Variant'], axis=1)
        df_final = df_final[cols_final]
        caller_variants_igv = pd.concat([caller_variants_igv,df_final])
    return caller_variants_igv

def somatic_annotation_filtered(all_caller_df,annotations_processed):
    final_df = pd.concat(all_caller_df)
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
    all_variants_final_df = all_variants_final[(all_variants_final['AF'] >= 0.05) & (all_variants_final['DP'] > 50) & (all_variants_final['ExonicFunc.refGene'] != 'synonymous SNV' ) & (all_variants_final['Func.refGene'] == 'exonic') | (all_variants_final['Func.refGene'] == 'upstream')| (all_variants_final['Func.refGene'] == 'splicing')]
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
    ## 1. copy files to igv external dirs ##
    copy_png_toexternal(rundir,runid)

    ## 2. process somatic variants ##
    somatic_annotation_file = find('annotations.paired.tsv', rundir)
    annotations_paired_file = pd.read_csv(somatic_annotation_file[0],sep="\t",header=0)
    annotations_processed = file_processing(annotations_paired_file)
    mutect2_variants_igv=strelka_variants_igv=lofreqsomatic_variants_igv=pd.DataFrame([])
    callers = ['MuTect2','Strelka','LoFreqSomatic']
    for caller in callers:
        variants_igv = add_igvsnapshots_to_variants(runid,rundir,caller,annotations_processed,"https://genome.med.nyu.edu/results/external/NYU/snuderllab/igv-snapshots/")
        if caller == "MuTect2":
            mutect2_variants_igv = variants_igv
        elif caller == "Strelka":
            strelka_variants_igv = variants_igv
        elif caller == "LoFreqSomatic":
            lofreqsomatic_variants_igv = variants_igv
        else:
            print("No caller found")

    all_caller_df = [mutect2_variants_igv,strelka_variants_igv,lofreqsomatic_variants_igv]
    somatic_variants = somatic_annotation_filtered(all_caller_df,annotations_processed)
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
    somatic_file_out = '%s/%s/%s'% (rundir,"clinical","somatic_variants.csv")
    somatic_variants_final.to_csv(somatic_file_out, index=False)

    somatic_variants_samples = somatic_variants_final.Test_Number.unique()

    ## 3. Process germline variants ##
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

    html_out_somatic = template_somatic.render(somatic_variants=somatic_variants_final,somatic_variants_samples=somatic_variants_samples,pactid=pactid,runid=runid)
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
