#!/usr/bin/env python3
import pandas as pd
import os
import argparse
import fnmatch
import numpy as np
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

def get_PACT_ID(demux_samplesheet_path):
    """get the previous run's PACT ID and run ID"""
    # Read the demux file
    try:
        with open(demux_samplesheet_path, 'r') as file:
            lines = file.readlines()
        # Find the row index where "[Data]" appears
        start_index = None
        for i, line in enumerate(lines):
            if line.startswith("[Data]"):
                start_index = i + 1
                break
        # "[Data]" row is found, read the CSV from the next row
        if start_index is not None:
            demux_data = pd.read_csv(demux_samplesheet_path, skiprows=start_index)
        else:
            print("No '[Data]' row found in the file.")
    except FileNotFoundError:
        print(f"File {demux_samplesheet_path} not found.")
        return
       
    PACT_ID = demux_data['Sample_Project'].unique()[0]
    Run_Number = demux_data['Run_Number'].unique()[0]
    return Run_Number, PACT_ID

def transform_filename(sample, runID, PACT_ID):
    """ Transform the sample filename for the boxplot """
    return sample.replace(runID, PACT_ID)

def determine_result(row):
    """ Determine the result based on the detected AF """
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
    # Step 3: Read the truth set data with AF range
    truth_set = pd.read_csv(sc_truthset, sep='\t')
    truth_set = truth_set.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    # Step 4: Merge SC truthset and current run sc annotation data
    combine_SCtruth_SCannot = pd.DataFrame()
    # For multiple SCs, handle the samples one by one to account for missing variants in only one sample
    for sampleID in sc_annot_data['Sample'].unique():
    # Merge tables and determine results
        sample_table = pd.merge(truth_set, sc_annot_data[sc_annot_data['Sample'] == sampleID], on=['Gene', 'Coding', 'AAChange'], how='left')
        sample_table['Sample'] = sampleID
        combine_SCtruth_SCannot=pd.concat([combine_SCtruth_SCannot, sample_table])
    combine_SCtruth_SCannot['Detected_AF'].fillna(0,inplace=True) # replace NA with 0
    combine_SCtruth_SCannot['Results'] = combine_SCtruth_SCannot.apply(determine_result, axis=1)
    return(combine_SCtruth_SCannot)

def get_directories(path):
    """Get a list of directory names in the given path."""
    directories = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    directories.sort()
    return directories

def get_previous_directories(directories, target_dir):
    """Get the target directory and the five preceding directories."""
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

def get_past5_results(runid, dir_NGS607,sc_vaf_table,demuxss):
    target_dir = runid
    previous_run_dirs = get_directories(dir_NGS607)
    previous_dir_results = get_previous_directories(previous_run_dirs, target_dir)
    print(previous_dir_results)
    # Merge data from previous runs with current run #
    if previous_dir_results:
        merged_prev_runs = sc_vaf_table[['Gene', 'Coding', 'AAChange', 'Detected_AF', 'Sample']]
        current_runInfo = get_PACT_ID(os.path.join(dir_NGS607, runid, demuxss))
        # get the pact id #
        current_PACTID = current_runInfo[1]
        print(current_PACTID) 
        # replace the samples with PACT ID of that run #
        merged_prev_runs['Sample'] = merged_prev_runs['Sample'].apply(lambda sample: transform_filename(sample, runid, current_PACTID))

        for prev_dir in previous_dir_results:
            file_path = os.path.join(dir_NGS607, prev_dir, "output/annotations/annotations.SeraCare.tsv")
            try:
                prev_data = pd.read_csv(file_path, sep='\t')
            except FileNotFoundError:
                 print(f"File {file_path} not found.")
                 continue
            prev_runInfo = get_PACT_ID(os.path.join(dir_NGS607, prev_dir, demuxss))
            prev_runID = prev_runInfo[0]
            prev_PACTID = prev_runInfo[1]
            prev_data = prev_data[['SeraCare.Gene', 'SeraCare.Coding', 'SeraCare.AAChange', 'AF', 'Sample']]
            prev_data.columns = ['Gene', 'Coding', 'AAChange', 'Detected_AF', 'Sample']
            prev_data['Sample'] = prev_data['Sample'].apply(lambda sample: transform_filename(sample, prev_runID, prev_PACTID))
            merged_prev_runs = pd.concat([merged_prev_runs, prev_data])
        else:
            print(f"Directory '{runid}' not found in the specified path.")
    # Apply the transformation and save the merged data
    return merged_prev_runs, current_PACTID


def plot_boxplot(previous_run_results,sc_output_past5runs_plot,current_PACTID):
    data = previous_run_results
    #Replace NA with 0
    data.fillna(0, inplace=True)
    # Extract part after the prefix "^[0-9]+_"
    data['extracted'] = data['Sample'].str.replace(r'^[0-9]+_', '', regex=True)
    # Set factor levels based on extracted strings
    sorted_samples = data.sort_values('extracted')['Sample'].unique()
    data['Sample'] = pd.Categorical(data['Sample'], categories=sorted_samples, ordered=True)
    # Assign colors based on pactid matching
    data['color'] = np.where(data['Sample'].str.contains(current_PACTID), 'lightgreen', 'coral')
    print(data)
    # Create the plot
    plot = (
    ggplot(data, aes(x='Sample', y='Detected_AF', fill='color')) +
    geom_boxplot(alpha=0.3) +
    stat_summary(aes(group='Sample'),fun_y=pd.Series.mean,geom='point',shape='o',size=2,color="red",fill="red")+
    theme(
        axis_text_x=element_text(angle=90, size=8, weight='bold'),
        axis_text_y=element_text(size=8, weight='bold'),
        plot_title=element_text(size=12, weight='bold'),
        axis_title_x=element_text(size=10, weight='bold'),
        axis_title_y=element_text(size=10, weight='bold'),
        legend_position='none'
    ) +
    labs(
        title='Distribution of SERACARE variants detected AF in current and past 5 runs',
        x='Sample',
        y='Detected AF'
    )
)
    # Save the plot
    output_path = sc_output_past5runs_plot
    plot.save(output_path,width=8,height=10,dpi=300)

def main(runid,rundir,sc_truthset,dir_NGS607):    
    # Get path for Seracare annotaiton file #
    seracare_annotation_file = find('annotations.SeraCare.tsv',rundir)
    # Get the output table for current run
    sc_vaf_table = process_run(seracare_annotation_file, sc_truthset)
    sc_output_file = '%s/%s/%s'% (rundir,"clinical","seracare_AF_QC_table.csv")
    sc_vaf_table.to_csv(sc_output_file,index=False)
    # # Get the output for past 5 runs
    demuxss ='demux-samplesheet.csv'
    previous_run_results, current_PACTID = get_past5_results(runid,dir_NGS607,sc_vaf_table,demuxss)
    # # Change the sample names for boxplot, write output for past 5 runs
    sc_output_past5runs_file = '%s/%s/%s'% (rundir,"clinical","past_5_runs_AF.csv")
    previous_run_results.to_csv(sc_output_past5runs_file,index=False)
    # Boxplot for VAF distribution across current and past 5 runs
    sc_output_past5runs_plot = '%s/%s/%s'% (rundir,"clinical","sc_AF_past_5_runs.png")
    plot_boxplot(previous_run_results,sc_output_past5runs_plot,current_PACTID)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-rid','--runid', help='Path to run id',required=True)
    parser.add_argument('-rdir','--rundir', help='Path to run directory',required=True)
    parser.add_argument('-sct','--sctruthset', help='SeraCare Ref file',required=True)
    parser.add_argument('-ngsdir','--ngs607dir',help='Path to all runs',required=True)
    args = parser.parse_args()
    main(args.runid,args.rundir,args.sctruthset,args.ngs607dir)