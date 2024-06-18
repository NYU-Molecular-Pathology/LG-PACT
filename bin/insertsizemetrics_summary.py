#!/usr/bin/env python3

import glob, os, argparse
import linecache
import pandas as pd
import jinja2

def insertmetrics_summary(rundir_path):
    # Get list of all files insertmetric txt #
    insertmetric_path_pattern = rundir_path+"output/CollectInsertSizeMetrics/*_insert_size_metrics.txt"
    insertmetrics_files_list = [f for f in glob.glob(insertmetric_path_pattern)]

    ##loop through each file, get sample name, line 7,8 ##
    all_dfs = []
    for file in insertmetrics_files_list:
        sampleID = file.split("/")[-1].split("_insert_size_metrics.txt")[0].split("_")[5]
        col_line = linecache.getline(file,7)
        col_line_lst = col_line.strip('\n').split('\t')
        col_line_lst.append("SampleID")
        val_line = linecache.getline(file,8)
        val_line_lst = val_line.strip('\n').split('\t')
        val_line_lst.append(sampleID)
        res = {col_line_lst[i]: val_line_lst[i] for i in range(len(col_line_lst))}
        df = pd.DataFrame(res,index=[0])
        all_dfs.append(df)
    merged_insertsizemetrics_df = pd.concat(all_dfs)
    return merged_insertsizemetrics_df

def insertmetrics_report(rundir_path):
    df_final = insertmetrics_summary(rundir_path)
    df_final_select_cols_for_htmlreport = df_final[['SampleID','MEDIAN_INSERT_SIZE','MODE_INSERT_SIZE','MEDIAN_ABSOLUTE_DEVIATION','MIN_INSERT_SIZE','MAX_INSERT_SIZE','MEAN_INSERT_SIZE','STANDARD_DEVIATION','READ_PAIRS']]
    insertmetrics_html_df = df_final_select_cols_for_htmlreport.sort_values(by=['SampleID'])

    template_path = os.path.join(rundir_path,"template/")
    templateLoader = jinja2.FileSystemLoader(searchpath=template_path)
    templateEnv = jinja2.Environment(loader=templateLoader)
    TEMPLATE_FILE = "insertmetricsreport_template.html"
    template_insertmetric = templateEnv.get_template(TEMPLATE_FILE)
    html_out_insertmetric = template_insertmetric.render(datains=insertmetrics_html_df)
    insertmetric_filename = '%s/%s'% (rundir_path,"output/clinical/insertmetrics-report.html")
    html_file_insertmetric = open(insertmetric_filename, 'w')
    html_file_insertmetric.write(html_out_insertmetric)
    html_file_insertmetric.close()

def main(rundir):
    insertmetrics_report(rundir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Provide path to run directory")
    parser.add_argument('-rdir', '--rundir', help='Your run dir path', required=True)
    args = parser.parse_args()
    main(args.rundir)