#!/usr/bin/env python3
import pandas as pd
import os
import argparse
import fnmatch
from plotnine import *
from plotnine.data import *
import patchworklib as pw

'''
Inputs:
inputfile - "annotations.seracare.tsv"
dir_NGS607 - path to NGS607 dir
sc-truthset - "SeraSeq-Truth.txt"
runID - RunID
rundir - Run output dir path

Outputs:
seracare_AF_QC_table.csv - VAF table for current run
seracare_AF_past_5_runs.csv - VAF table for past 5 runs
sc_AF_past_5_runs.png - Boxplot VAF distribution for current and past 5 runs

'''

def find(pattern, path):
    exclude_directories = set(['work'])
    file_path = []
    for root,dirs,files in os.walk(path):
        dirs[:] = [d for d in dirs if d not in exclude_directories]
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                file_path.append(os.path.join(root, name))
    return file_path

def determine_result(row):
    lower_limit = float(row['AF 95% Interval (Historical Runs)'].split('~')[0].strip())
    if row['Detected_AF'] == 0:
        return "Not detected"
    elif row['Detected_AF'] <= lower_limit:
        return "Low"
    else:
        return "Yes"
    
def process_run(seracare_annotation_file, sc_truthset):
    # Step 1: Read the TSV file
    sc_annot_data = pd.read_csv(seracare_annotation_file[0], sep='\t')
    # Step 2: Select the necessary columns
    sc_annot_data = sc_annot_data[['SeraCare.Gene', 'SeraCare.Coding', 'SeraCare.AAChange', 'DP', 'AF', 'Sample']]
    sc_annot_data.columns = ['Gene', 'Coding', 'AAChange', 'DP', 'Detected_AF', 'Sample'] # include samples to account for multiple seracares per run
    # Step 3: Read the truth set data
    truth_set = pd.read_csv(sc_truthset, sep='\t')
    truth_set = truth_set.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    # Step 4: Merge SC truthset and current run sc annotation data
    combine_SCtruth_SCannot = pd.merge(truth_set, sc_annot_data, on=['Gene', 'Coding', 'AAChange'], how='left')
    # If 'Sample' column has NaN, fill it with the sample value from 'data'
    if combine_SCtruth_SCannot['Sample'].isna().any():
        if len(sc_annot_data['Sample'].unique()) == 1:
            combine_SCtruth_SCannot['Sample'].fillna(sc_annot_data['Sample'].iloc[0], inplace=True)
        else:
            combine_SCtruth_SCannot['Sample'] = combine_SCtruth_SCannot.groupby(['Gene', 'Coding', 'AAChange'])['Sample'].transform(lambda x: x.ffill().bfill())
    print(combine_SCtruth_SCannot.to_string())
    combine_SCtruth_SCannot['Detected_AF'] = combine_SCtruth_SCannot['Detected_AF'].fillna(value=0) # replace NA with 0
    combine_SCtruth_SCannot['Results'] = combine_SCtruth_SCannot.apply(determine_result, axis=1)
    return(combine_SCtruth_SCannot)

def get_directories(path):
    """Get a list of directory names in the given path."""
    directories = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    directories.sort()
    return directories

def get_previous_directories(directories, target_dir):
    """Get the target directory and the four preceding directories."""
    if target_dir in directories:
        index = directories.index(target_dir)
        preceding_dirs = []
        required_count = 5 # required 5 runs
        seen_pact = set()

        # Iterate backwards from the target run's index
        for i in range(index - 1, -1, -1):
            dir = directories[i]
            parts = dir.split('_')
            if len(parts) > 4:  # Treat as having extra suffix
                base_dir = '_'.join(parts[:4])
            else:
                base_dir = dir

            if base_dir not in seen_pact:
                seen_pact.add(base_dir)
                preceding_dirs.insert(0, dir)
                required_count -= 1
                if required_count == 0:
                    break

        return preceding_dirs
    else:
        return None

def get_past5_results(runid, dir_NGS607,sc_vaf_table):
    target_dir = runid
    print(target_dir)
    previous_run_dirs = get_directories(dir_NGS607)
    previous_dir_results = get_previous_directories(previous_run_dirs, target_dir)
    print(previous_dir_results)
    if previous_dir_results:
        all_runs_merged = sc_vaf_table[['Gene', 'Coding', 'AAChange', 'Detected_AF', 'Sample']]
        for result in previous_dir_results:
            file_path = os.path.join(dir_NGS607, result, "output/annotations/annotations.SeraCare.tsv")
            try:
                data = pd.read_csv(file_path, sep='\t')
            except FileNotFoundError:
                print(f"File {file_path} not found.")
            previour_dirs_data = data[['SeraCare.Gene', 'SeraCare.Coding', 'SeraCare.AAChange', 'AF', 'Sample']]
            previour_dirs_data.columns = ['Gene', 'Coding', 'AAChange', 'Detected_AF', 'Sample']
            all_runs_merged = pd.concat([all_runs_merged, previour_dirs_data])
        return all_runs_merged
    else:
        print(f"Directory '{target_dir}' not found in the specified path.")

def transform_filename(sample):
    parts = sample.split('_')
    new_samplname = f"{parts[1]}_{parts[2]}_{parts[3]}_{parts[4]}_{parts[0]}"
    return new_samplname

def plot_boxplot(previous_run_results,sc_output_past5runs_plot):
    data = previous_run_results
    plot = (ggplot(data, aes(x="Sample", y="Detected_AF", fill="Sample"))+geom_boxplot(alpha=0.3)+scale_fill_manual(values=["red","red","red","red","red","green"])+theme_bw()+theme(axis_text_x=element_text(angle=90,hjust=1))+labs(title="Distribution of SERACARE variants detected AF in current and past 5 runs"))
    boxplot_past5_runs = pw.load_ggplot(plot, figsize=(15,10))
    return boxplot_past5_runs.savefig(sc_output_past5runs_plot)

def main(runid,rundir,sc_truthset,dir_NGS607):    
    # Get path for Seracare annotaiton file #
    seracare_annotation_file = find('annotations.SeraCare.tsv',rundir)
    # Get the output table for current run
    sc_vaf_table = process_run(seracare_annotation_file, sc_truthset)
    sc_output_file = '%s/%s/%s'% (rundir,"clinical","seracare_AF_QC_table.csv")
    sc_vaf_table.to_csv(sc_output_file,index=False)
    # Get the output for past 5 runs
    previous_run_results = get_past5_results(runid,dir_NGS607,sc_vaf_table)
    # Change the sample names for boxplot, write output for past 5 runs
    previous_run_results['Sample'] = previous_run_results['Sample'].apply(transform_filename)
    sc_output_past5runs_file = '%s/%s/%s'% (rundir,"clinical","seracare_AF_past_5_runs.csv")
    previous_run_results.to_csv(sc_output_past5runs_file,index=False)
    # Boxplot for VAF distribution across current and past 5 runs
    sc_output_past5runs_plot = '%s/%s/%s'% (rundir,"clinical","sc_AF_past_5_runs.png")
    plot_boxplot(previous_run_results,sc_output_past5runs_plot)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-rid','--runid', help='Path to run id',required=True)
    parser.add_argument('-rdir','--rundir', help='Path to run directory',required=True)
    parser.add_argument('-sct','--sctruthset', help='SeraCare Ref file',required=True)
    parser.add_argument('-ngsdir','--ngs607dir',help='Path to all runs',required=True)
    args = parser.parse_args()
    main(args.runid,args.rundir,args.sctruthset,args.ngs607dir)