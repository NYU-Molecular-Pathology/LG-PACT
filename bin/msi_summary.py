#!/usr/bin/env python3
"""
Generate MSI summary from msisensor results.

"""

import argparse
import pandas as pd

def read_samplesheet(samplesheet_path):
    """Read and filter the samplesheet to get the T:N filename pairs list."""
    samplesheet = pd.read_csv(samplesheet_path, skiprows=19)
    samplesheet_paired = samplesheet.dropna()
    samplesheet_filtered = samplesheet_paired[
        ~samplesheet_paired['Sample_Name'].str.contains("SC_SERACARE|NC_HAPMAP|NTC_H20")
    ]
    samplesheet_filtered['Sample_pairs'] = samplesheet_filtered["Sample_Name"] + "_" + samplesheet_filtered["Paired_Normal"]
    return samplesheet_filtered

def get_msi_data(sample_list, run_dir):
    """Retrieve MSI data for all sample pairs from the MSI results directory."""
    all_data = []
    for sample in sample_list:
        msi_path = f"{run_dir}/msi/{sample}"
        msi_data = pd.read_csv(msi_path, sep="\t")
        msi_data['Sample_pairs'] = sample
        all_data.append(msi_data)
    return pd.concat(all_data)

def merge_dicts_to_summary(demux_ss_dict, msi_dict):
    """Merge dictionaries to create the final MSI summary."""
    return {k: tuple(d[k] for d in [demux_ss_dict, msi_dict]) for k in msi_dict.keys()}

def write_msi_summary(msi_summary, output_file):
    """Write the MSI summary to an output file."""
    with open(output_file, 'w') as fout:
        fout.write("SampleID\tTotal_Number_of_Sites\tNumber_of_SomaticSites\tMSI\tStatus\tTumorType\n")
        for samplepairs, msi_results in msi_summary.items():
            SampleID = "_".join(samplepairs.split("_")[:7])
            TumorType = msi_results[0]
            Total_Sites, Somatic_sites = map(int, msi_results[1][:2])
            msi = "NA" if Total_Sites == 0 else msi_results[1][2]
            msi_status = "NA" if msi == "NA" else ("Unstable" if msi >= 14.0 else "Stable")
            fout.write(f"{SampleID}\t{Total_Sites}\t{Somatic_sites}\t{msi}\t{msi_status}\t{TumorType}\n")

def get_options():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-rdir', '--rundir', help='Your run dir path', required=True)
    parser.add_argument("-s", "--sample_sheet", type=str, required=True, help="Demux sample sheet")
    return parser.parse_args()

def main():
    """Main function to generate the MSI summary."""
    args = get_options()
    demux_ss = read_samplesheet(args.sample_sheet)
    sample_list = demux_ss['Sample_pairs'].tolist()
    demux_ss_dict = dict(zip(demux_ss.Sample_pairs, demux_ss.Tumor_Type))
    
    msi_df = get_msi_data(sample_list, args.rundir)
    msi_df.rename(columns={"%": "MSI"}, inplace=True)
    msi_dict = msi_df.set_index('Sample_pairs').T.to_dict('list')
    
    msi_summary = merge_dicts_to_summary(demux_ss_dict, msi_dict)
    msi_summary_file = '%s/%s/%s'% (args.rundir,"clinical","msi_validation.tsv")
    write_msi_summary(msi_summary, msi_summary_file)

if __name__ == "__main__":
    main()
