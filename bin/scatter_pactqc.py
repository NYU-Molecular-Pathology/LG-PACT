#!/usr/bin/env python3

import pandas as pd
import os,fnmatch,seaborn
from functools import reduce
import math, argparse
from plotnine import *
from plotnine.data import *
import patchworklib as pw
import matplotlib.cm
import matplotlib
from matplotlib import cm as colormaps

### Find annot file for all 3 callers ###
def find(pattern, path):
    exclude_directories = set(['work','variants','alignments'])
    file_path = []
    for root,dirs,files in os.walk(path):
        dirs[:] = [d for d in dirs if d not in exclude_directories]
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                file_path.append(os.path.join(root, name))
    return file_path


def annots_3callers(file):
    data = pd.read_csv(file, sep="\t")
    caller = file.split('/')[-1].split(".")[1]
    data_filter = data.loc[(data["AF"] >= 0.04)]
    data_variants = data_filter.copy()
    data_variants["VariantFinal"] = data_variants['CHROM']+':'+ data_variants['POS'].astype(str)+':'+data_variants['REF']+':'+data_variants['ALT']+':'+data_variants['Tumor']
    reqd_cols = data_variants[['VariantFinal','QUAL','AF','DP','Tumor','Normal','Func.refGene','ExonicFunc.refGene','Gene.refGene','AAChange.refGene']]
    suffix = "_"+caller
    reqd_cols.columns = reqd_cols.columns.map(lambda x : x+suffix if x != 'VariantFinal' else x)
    return reqd_cols


def plot_qc(somatic_allcaller_results,pactid,rundir):
    somatic_allcaller_results['QUAL_LoFreqSomatic'] = somatic_allcaller_results['QUAL_LoFreqSomatic'].apply(lambda x: math.sqrt(float(x)))
    p1_ggtitle = pactid+".A.3callers"
    p1 = (ggplot(somatic_allcaller_results,aes(y="QUAL_Strelka", x="QUAL_MuTect2",size="QUAL_LoFreqSomatic",color='QUAL_LoFreqSomatic'))+geom_point()+facet_wrap('~Test_Number')+ggtitle(p1_ggtitle)+theme_bw()+labs(color="sqrt(QUAL.LFS)",size="sqrt(QUAL.LFS)"))
    g1 = pw.load_ggplot(p1, figsize=(15,10))

    somatic_allcaller_results['QUAL_Strelka'] = somatic_allcaller_results['QUAL_Strelka'].apply(lambda x: math.log2(x))
    somatic_allcaller_results['QUAL_MuTect2'] = somatic_allcaller_results['QUAL_MuTect2'].apply(lambda x: math.log2(x))
    p2_ggtitle = pactid+".B.3callers"
    p2 = (ggplot(somatic_allcaller_results,aes(y="QUAL_Strelka", x="QUAL_MuTect2",size="QUAL_LoFreqSomatic",color='QUAL_LoFreqSomatic'))+geom_point()+facet_wrap('~Test_Number')+ggtitle(p2_ggtitle)+theme_bw()+labs(x="log2(QUAL_MuTect2)",y="log2(QUAL_Strelka)",color="sqrt(QUAL.LFS)",size="sqrt(QUAL.LFS)"))
    g2 = pw.load_ggplot(p2, figsize=(15,10))
    g12 = (g1/g2)
    figure_title = "PACT QC "+ pactid
    output_path = rundir+"/clinical/"+figure_title
    return g12.savefig(output_path)


def main(**kwargs):
    rundir = kwargs.pop('rundir', None)
    pactid = kwargs.pop('pactid', None)

    ### process annotation files for Mutect2, Strelka and Lfs ###
    files = ["annotations.MuTect2.tsv","annotations.Strelka.tsv","annotations.LoFreqSomatic.tsv"]
    annotations_dir = rundir+"/annotations/"
    print("Processing annotation files!!")
    annotation_file = [find(file, annotations_dir) for file in files]
    somatic_variants_callers = [annots_3callers(caller[0]) for caller in annotation_file ]
    final_df = reduce(lambda  left,right: pd.merge(left,right,on='VariantFinal',how='outer'), somatic_variants_callers).fillna('NA')
    final_df['Sample'] = final_df['VariantFinal'].str.split(":",n = 4, expand = True)[4]
    final_df[['QUAL_MuTect2', 'QUAL_Strelka','QUAL_LoFreqSomatic']] = final_df[['QUAL_MuTect2','QUAL_Strelka','QUAL_LoFreqSomatic']].replace('NA',value=0.81)

    ## merge the ngs test case ID with the processed annotaiton of variants##
    maindir = rundir.rpartition('/')[0]
    demux_ss_find = find('demux-samplesheet.csv', maindir)
    demux_ss_file = pd.read_csv(demux_ss_find[0],skiprows=19)
    somatic_demux_merge = pd.merge(final_df,demux_ss_file,left_on="Sample",right_on="Sample_ID")
    somatic_demux_merge_excludedHAPMAP = somatic_demux_merge.loc[(somatic_demux_merge['Test_Number'] != '0')]
    variants_final = somatic_demux_merge_excludedHAPMAP[["Test_Number","VariantFinal" ,"QUAL_MuTect2", "AF_MuTect2","DP_MuTect2","Tumor_MuTect2",   "Normal_MuTect2","Func.refGene_MuTect2","ExonicFunc.refGene_MuTect2", "Gene.refGene_MuTect2","AAChange.refGene_MuTect2","QUAL_Strelka","AF_Strelka",  "DP_Strelka","Tumor_Strelka","Normal_Strelka","Func.refGene_Strelka","ExonicFunc.refGene_Strelka","Gene.refGene_Strelka","AAChange.refGene_Strelka","QUAL_LoFreqSomatic","AF_LoFreqSomatic","DP_LoFreqSomatic","Tumor_LoFreqSomatic","Normal_LoFreqSomatic","Func.refGene_LoFreqSomatic","ExonicFunc.refGene_LoFreqSomatic","Gene.refGene_LoFreqSomatic","AAChange.refGene_LoFreqSomatic","Sample","Sample_ID","Sample_Name","Paired_Normal","I7_Index_ID","index","Specimen_ID","EPIC_ID","Tumor_Content","Tumor_Type","Description","Run_Number","Sequencer_ID","Chip_ID","Sample_Project"]]
    output_path = rundir+"/clinical/"+"PACT.3callers.csv"
    variants_final.to_csv(output_path, index=False)
    Pact_3callers = variants_final[['Test_Number','VariantFinal','QUAL_MuTect2','QUAL_Strelka','QUAL_LoFreqSomatic']]

    ## plot the qual of 3 callers (based of xf) ##
    print("Generating QUAL figure and PACT 3 callers csv file!!")
    plot_qc(Pact_3callers,pactid,rundir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Provide path to run directory's output dir")
    parser.add_argument('-rdir','--rundir', help='Your run directory',required=True)
    parser.add_argument('-pactid','--pactid',help='PACT ID',required=True )
    args = parser.parse_args()
    main(**vars(args))
