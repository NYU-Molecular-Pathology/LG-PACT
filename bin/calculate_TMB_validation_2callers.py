#!/usr/bin/env python3
"""
Calculate the tumor mutation burden in variants/Megabase for Validation (MuTect2+LoFreqSomatic)

"""
import argparse
import pandas as pd

"""
Additional criteria:

For both matched and unmatched we will apply the following criteria:
1- VAF >5%
2- coverage > 200X.
3- include Non Synonymous and synonymous SNV.
4- include only exonic.
5- variants called by both MuTect and LoFreqSomatic only

"""
def load_data(demux_ss, annot_paired):
    demux_ss_df = pd.read_csv(demux_ss, skiprows=19)
    annot_paired_df = pd.read_csv(annot_paired, sep="\t")
    return demux_ss_df, annot_paired_df

def filter_annotation_data(annot_paired_df):
    annot_paired_data = annot_paired_df.astype({'POS': str})
    annot_paired_filter_AFDP = annot_paired_data[(annot_paired_data['AF'] > 0.05) & (annot_paired_data['DP'] >= 200)]
    annot_paired_type_filter = annot_paired_filter_AFDP[((annot_paired_filter_AFDP['Func.refGene'] == "exonic") & (annot_paired_filter_AFDP['ExonicFunc.refGene'].isin(["nonsynonymous SNV", "synonymous SNV"])))]
    annot_paired_type_filter['Variant'] = annot_paired_type_filter[['CHROM', 'POS', 'REF', 'ALT']].apply(lambda x: ':'.join(x), axis=1)
    annot_Mutect_LoFreq = annot_paired_type_filter[(annot_paired_type_filter['VariantCaller'].isin(["MuTect2", "LoFreqSomatic"]))]
    annot_Mutect_LoFreq_nocontrols = annot_Mutect_LoFreq[~annot_Mutect_LoFreq['Tumor'].str.contains("SC-SERACARE|NC-HAPMAP|SC_SERACARE|NC_HAPMAP")]
    return annot_Mutect_LoFreq_nocontrols

def filter_demux_data(demux_ss_df):
    demux_ss_data = demux_ss_df.dropna()
    demux_ss_data_final = demux_ss_data[~demux_ss_data['Sample_Name'].str.contains("SC_SERACARE|NC_HAPMAP|NTC_H20")]
    demux_ss_dict = dict(zip(demux_ss_data_final.Sample_Name, demux_ss_data_final.Tumor_Type))
    return demux_ss_data_final, demux_ss_dict

def compute_variant_counts(samples, annot_Mutect_LoFreq_nocontrols):
    sample_nvariants_dict = {}
    for each_sample in samples:
        sample_df = annot_Mutect_LoFreq_nocontrols[annot_Mutect_LoFreq_nocontrols['Tumor'] == each_sample]
        mutect_variants = sample_df[sample_df['VariantCaller'] == "MuTect2"]
        lofreqsomatic_variants = sample_df[sample_df['VariantCaller'] == "LoFreqSomatic"]
        common_variants = pd.merge(mutect_variants, lofreqsomatic_variants, how='inner', on=['Variant'])
        n_common_variants = len(common_variants)
        sample_nvariants_dict[each_sample] = n_common_variants
    return sample_nvariants_dict

def load_callable_loci(callable_loci):
    loci_dict = {}
    with open(callable_loci) as fl:
        for line in fl:
            items = line.strip().split()
            loci_dict[items[1]] = items[0]
    return loci_dict

def merge_dicts(sample_nvariants_dict, loci_dict, demux_ss_dict):
    dict_list = [sample_nvariants_dict, loci_dict, demux_ss_dict]
    tmb_vals_dict = {k: tuple(d[k] for d in dict_list) for k in sample_nvariants_dict.keys()}
    return tmb_vals_dict

def write_output(tmb_vals_dict, output_tmb_file):
    with open(output_tmb_file, 'w') as fout:
        fout.write("SampleID\tVariantCaller\tnBases\tnVariants\tTMB\tTumorType\n")
        for sampleID, values in tmb_vals_dict.items():
            if values[1] == "0":
                tmb = "NA"
                fout.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (sampleID, "MuTect2+LoFreqsomatic", str(values[1]),
                                                     str(values[0]), str(tmb), str(values[2])))
            else:
                tmb = round(float(values[0]) / float(values[1]) * 1000000, 2)
                fout.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (sampleID, "MuTect2+LoFreqsomatic", str(values[1]),
                                                     str(values[0]), str(tmb), str(values[2])))
               
def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--loci", type=str, required=True,
                        help="callable loci from alignments")
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Path to input variant annotation including file name")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Path to output TMB file including file name")
    parser.add_argument("-s", "--sample_sheet", type=str, required=True, help="sample sheet")
    parser.add_argument("-t", "--type", type=str, required=False,
                        help="tumor sample type, paired not unpaired", default="paired")
    return parser.parse_args()

def main():
    args = get_options()
    # Load demux ss and annotations paired data
    demux_ss_df, annot_paired_df = load_data(args.sample_sheet, args.input)
    # Apply filter criteria for TMB
    annot_Mutect_LoFreq_nocontrols = filter_annotation_data(annot_paired_df)
    # Get sample and tumor type list
    demux_ss_data_final, demux_ss_dict = filter_demux_data(demux_ss_df)
    samples = demux_ss_data_final.Sample_Name.unique()
    # Get variant count for each sample after applying filters
    sample_nvariants_dict = compute_variant_counts(samples, annot_Mutect_LoFreq_nocontrols)
    # Load the callable loci data
    loci_dict = load_callable_loci(args.loci)
    # Calculate TMB
    tmb_vals_dict = merge_dicts(sample_nvariants_dict, loci_dict, demux_ss_dict)
    # Write final TMB output
    write_output(tmb_vals_dict,args.output)

if __name__ == "__main__":
    main()