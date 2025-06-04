#!/usr/bin/env python3

import argparse
import pandas as pd
import glob

def parse_contamination(rundir_path):
    # Get list of all contamination.txt #
    conpair_contamination_path_pattern = rundir_path+"output/Conpair/*_contamination.txt"
    contamination_files_list = [f for f in glob.glob(conpair_contamination_path_pattern)]
    contamination_df = []
    for contamination_file in contamination_files_list:
        sampleID = contamination_file.split("/")[-1].split("_contamination.txt")[0]
        with open(contamination_file) as f:
            lines = f.readlines()
            normal = None
            tumor = None
            for line in lines:
                if "Normal sample contamination level" in line:
                    val_normal = float(line.split(":")[1].strip().replace('%', ''))
                    normal = "<0.1" if val_normal < 0.1 else val_normal
                elif "Tumor sample contamination level" in line:
                    val_tumor = float(line.split(":")[1].strip().replace('%', ''))
                    tumor = "<0.1" if val_tumor < 0.1 else val_tumor
        contamination_df.append([sampleID, normal, tumor])
    return pd.DataFrame(contamination_df, columns=["SampleID", "Normal_contamination", "Tumor_contamination"])

def parse_concordance(rundir_path):
    # Get list of all concordance.txt #
    conpair_concordance_path_pattern = rundir_path+"output/Conpair/*_concordance.txt"
    concordance_files_list = [f for f in glob.glob(conpair_concordance_path_pattern)]
    concordance_df = []
    for concordance_file in concordance_files_list:
        sampleID = concordance_file.split("/")[-1].split("_concordance.txt")[0]
        with open(concordance_file) as f:
            lines = f.readlines()
            concordance = None
            for line in lines:
                if "Concordance" in line:
                    concordance = float(line.split(":")[1].strip().replace('%', ''))
        concordance_df.append([sampleID, concordance])
    return pd.DataFrame(concordance_df, columns=["SampleID", "Conpair_Concordance"])

def combine_conpair_data(rundir_path):
    # combine the summary
    contamination_df = parse_contamination(rundir_path)
    concordance_df = parse_concordance(rundir_path)
    conpair_summary = pd.merge(contamination_df, concordance_df, on="SampleID", how="outer")
    conpair_summary_noHapMap = conpair_summary[~conpair_summary["SampleID"].str.contains("HapMap", na=False)]
    conpair_output_path = rundir_path+"output/clinical/"+"conpair_summary.csv"
    conpair_summary_noHapMap.to_csv(conpair_output_path,index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Provide path to run directory")
    parser.add_argument('-rdir', '--rundir', help='Your run dir path', required=True)
    args = parser.parse_args()
    combine_conpair_data(args.rundir)
