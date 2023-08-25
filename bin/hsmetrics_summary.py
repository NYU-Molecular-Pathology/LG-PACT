#!/usr/bin/env python3

import glob, os, argparse
import linecache
import pandas as pd
import jinja2
from plotnine import *
from plotnine.data import *
import patchworklib as pw

def hsmetrics_summary(rundir_path):
    # Get list of all files hsmetric txt #
    hsmetric_path_pattern = rundir_path+"output/CollectHsMetrics/*_hs_metrics.txt"
    hsmetric_files_list = [f for f in glob.glob(hsmetric_path_pattern)]
    ##loop through each file, get sample name, line 7,8 ##
    all_dfs = []
    for file in hsmetric_files_list:
        sampleID = file.split("/")[-1].split("_hs_metrics.txt")[0].split("_")[5]
        col_line = linecache.getline(file,7)
        col_line_lst = col_line.strip('\n').split('\t')
        col_line_lst.append("SampleID")
        val_line = linecache.getline(file,8)
        val_line_lst = val_line.strip('\n').split('\t')
        val_line_lst.append(sampleID)
        res = {col_line_lst[i]: val_line_lst[i] for i in range(len(col_line_lst))}
        df = pd.DataFrame(res,index=[0])
        all_dfs.append(df)
    merged_hsmetric_df = pd.concat(all_dfs)
    return merged_hsmetric_df


def hsmetric_report(rundir_path):
    df_final = hsmetrics_summary(rundir_path)
    df_final_select_cols_for_htmlreport = df_final[['SampleID','ON_BAIT_BASES','MEAN_BAIT_COVERAGE','PCT_USABLE_BASES_ON_BAIT','PCT_USABLE_BASES_ON_TARGET','TOTAL_READS','PF_READS','ON_TARGET_BASES','MEAN_TARGET_COVERAGE','FOLD_80_BASE_PENALTY','AT_DROPOUT','GC_DROPOUT']]
    hsmetric_html_df = df_final_select_cols_for_htmlreport.sort_values(by=['SampleID'])

    template_path = os.path.join(rundir_path,"template/")
    templateLoader = jinja2.FileSystemLoader(searchpath=template_path)
    templateEnv = jinja2.Environment(loader=templateLoader)
    TEMPLATE_FILE = "hsreport_template.html"
    template_hs = templateEnv.get_template(TEMPLATE_FILE)
    html_out_hsmetric = template_hs.render(datahs=hsmetric_html_df)
    hs_filename = '%s/%s'% (rundir_path,"output/clinical/hsmetrics-report.html")
    html_file_hsmetric = open(hs_filename, 'w')
    html_file_hsmetric.write(html_out_hsmetric)
    html_file_hsmetric.close()

def lineplot(rundir_path,tumor_normal_csv):
    df_final = hsmetrics_summary(rundir_path)
    df_final_select_cols_for_plot = df_final[["SampleID","PCT_TARGET_BASES_10X","PCT_TARGET_BASES_20X","PCT_TARGET_BASES_30X","PCT_TARGET_BASES_40X","PCT_TARGET_BASES_50X","PCT_TARGET_BASES_100X","PCT_TARGET_BASES_250X","PCT_TARGET_BASES_500X","PCT_TARGET_BASES_1000X"]]
    sample_tumor_normal = pd.read_csv(tumor_normal_csv)
    tumor_samples = sample_tumor_normal['Tumor'].str.split("_",expand=True)[5].astype(str)
    normal_samples = sample_tumor_normal['Normal'].str.split("_",expand=True)[5].astype(str)

    df_pct_normal = df_final_select_cols_for_plot[df_final_select_cols_for_plot.SampleID.isin(normal_samples)]
    df_pct_normal_sort = df_pct_normal.sort_values(by=['SampleID'])
    df_pct_normal_melt = pd.melt(df_pct_normal_sort, id_vars=['SampleID'])
    df_pct_normal_melt.columns = ["SampleID","Group","Value"]
    df_pct_normal_melt['Value'] = df_pct_normal_melt['Value'].astype(float).round(2)
    normal_plt = (ggplot(df_pct_normal_melt,aes(x="Group",y="Value",color="SampleID",group="SampleID"))+geom_point()+geom_line()+scale_x_discrete(limits=("PCT_TARGET_BASES_10X","PCT_TARGET_BASES_20X","PCT_TARGET_BASES_30X","PCT_TARGET_BASES_40X","PCT_TARGET_BASES_50X","PCT_TARGET_BASES_100X","PCT_TARGET_BASES_250X","PCT_TARGET_BASES_500X","PCT_TARGET_BASES_1000X"))+geom_vline(xintercept = 7, linetype="dotted", color = "red")+theme_classic())
    normal_plt_gp = pw.load_ggplot(normal_plt, figsize=(15,10))
    
    df_pct_tumor = df_final_select_cols_for_plot[df_final_select_cols_for_plot.SampleID.isin(tumor_samples)]
    df_pct_tumor_sort = df_pct_tumor.sort_values(by=['SampleID'])
    df_pct_tumor_melt = pd.melt(df_pct_tumor_sort, id_vars=['SampleID'])
    df_pct_tumor_melt.columns = ["SampleID","Group","Value"]
    df_pct_tumor_melt['Value'] = df_pct_tumor_melt['Value'].astype(float).round(2)
    tumor_plt = (ggplot(df_pct_tumor_melt,aes(x="Group",y="Value",color="SampleID",group="SampleID"))+geom_point()+geom_line()+scale_x_discrete(limits=("PCT_TARGET_BASES_10X","PCT_TARGET_BASES_20X","PCT_TARGET_BASES_30X","PCT_TARGET_BASES_40X","PCT_TARGET_BASES_50X","PCT_TARGET_BASES_100X","PCT_TARGET_BASES_250X","PCT_TARGET_BASES_500X","PCT_TARGET_BASES_1000X"))+geom_vline(xintercept = 7, linetype="dotted", color = "red")+theme_classic())
    tumor_plt_gp = pw.load_ggplot(tumor_plt, figsize=(15,10))
    g12 = (normal_plt_gp/tumor_plt_gp)
    figure_title = "hsmetric_lineplot.png"
    output_path = rundir_path+"output/clinical/"+figure_title
    return g12.savefig(output_path)

def hsmetric_targetcoverage(rundir_path):
    pertarget_path_pattern = rundir_path+"output/CollectHsMetrics/*PerTargetCoverage.txt"
    pertarget_files_list = [f for f in glob.glob(pertarget_path_pattern)]
    processed_dfs = []
    for targcov_file in pertarget_files_list:
        filename = targcov_file.split("/")[-1]
        samplename = filename.split("_PerTargetCoverage.txt")[0]
        samplename_final = samplename.split("_")[5]
        process_data = pd.read_csv(targcov_file,sep="\t")
        process_data['SampleID'] = samplename_final
        process_readcount_filter = process_data[(process_data['read_count'] < 200)]
        data_final_df = process_readcount_filter[['SampleID','chrom','start','end','name','mean_coverage','normalized_coverage','read_count']]
        processed_dfs.append(data_final_df)
    perTargCovg_probes = pd.concat(processed_dfs).sort_values(by=["SampleID"])
    perTargCovg_probes_noNTC = perTargCovg_probes[(perTargCovg_probes["SampleID"]!= "NTC")]
    pertarg_output_path = rundir_path+"output/clinical/"+"PerTargetCoverage.csv"
    perTargCovg_probes_noNTC.to_csv(pertarg_output_path,index=False)


def main(rundir,tumor_normal_csv):
    hsmetric_report(rundir)
    lineplot(rundir,tumor_normal_csv)
    hsmetric_targetcoverage(rundir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Provide path to run directory")
    parser.add_argument('-rdir', '--rundir', help='Your run dir path', required=True)
    parser.add_argument('-tnss', '--tnss', help='Path to samples.tumor.normal.csv file', required=True)
    args = parser.parse_args()
    main(args.rundir,args.tnss)