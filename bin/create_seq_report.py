#!/usr/bin/env python

__author__ = "Kelsey Zhu"
__version__ = "1.0.1"

import pandas as pd
import os

class seq_qc_report:
    """
    Total number of reads > 10 million.
    Number of deduplicated mapped reads >10 million
    Number of targets with <50 coverage < 1%
    Average and median coverage of the sample >300X
    """

    def __init__(self, args):
        self.args = args

    # check total number of reads
    def check_total_num_of_reads(self, row, cutoff=10000000):
        mapping_reads = row['TotalMappedReads']
        if "_NTC_" in row["Sample"] or "_NTC2_" in row["Sample"]:
            return "%s Pass"%mapping_reads if mapping_reads < 30000 else "%s Fail"%mapping_reads
        else:
            return "%s Pass"%mapping_reads if mapping_reads > cutoff else "%s Fail"%mapping_reads

    # check sample coverage
    def check_sample_coverage(self, row, cutoff=0.01):
        pcnt = 100-float(row['50'])
        if pcnt < cutoff*100:
            return "%s"%str(round(pcnt,2))+"% Pass"
        else:
            return "%s"%str(round(pcnt,2))+"% Fail"

     # check duplication rate    
    def check_dup_rate(self, row, cutoff=0.50):
        dedup_rate = row['Deduplication_Rate']
        if float( dedup_rate < cutoff):
            return "%s"%str(round(dedup_rate,2))+" Pass"
        else:
            return "%s"%str(round(dedup_rate,2))+" Fail"

    # Average and median coverage of the sample
    def check_average_median_coverage(self, row, cutoff=300):
        sample = row["Sample"]
        if "_NTC_" in sample or "_NTC2_" in sample:
            return "%s Pass"%row["MedianCoverage"] if row["MedianCoverage"] < 5 else "%s Fail"%row["MedianCoverage"]
        elif sample.endswith("HAPMAP"):
            return "%s Pass"%row["MeanCoverage"] if row["MeanCoverage"] > cutoff else "%s Fail"%row["MeanCoverage"]
        elif row["MeanCoverage"] > cutoff and row["MedianCoverage"] > cutoff:
            return "%s,%s Pass" %(row["MeanCoverage"],row["MedianCoverage"])
        else:
            return "%s,%s Fail" %(row["MeanCoverage"],row["MedianCoverage"])

    def check_snp_overlap(self, row, cutoff=90):
        sample = row["Sample"]
        if "_NTC_" in sample or "_NTC2_" in sample or sample.endswith(tuple(excl_list)):
            return "NA"
        overlap = snp_overlap.loc[(snp_overlap['Tumor'] == sample) | (snp_overlap['Normal'] == sample)]
        try:
            pcnt = overlap.iloc[0]['SNP_Overlap_Pcent(%)']
            total_snps = overlap.iloc[0]['N_SNP_Count'].astype(str)
            return "%s"%pcnt+'%' + ' (out of %s) Pass'%total_snps if pcnt > cutoff else "%s"%pcnt+'%' + ' (out of %s) Fail'%total_snps
        except:
            return "NA"

    def check_negative_control_VAF(self, row, cutoff=0.02):
        sample = row["Sample"]
        if not sample.endswith("HAPMAP"):
            return "NA"
        var_rows = len(hapmap_variants.index)
        total_rows = len(unpaired_variants.loc[unpaired_variants['Sample'] == sample].index)
        vaf = round(float(var_rows/total_rows),2)
        return "%s"%str(vaf*100)+'% Pass' if vaf <= cutoff else "%s"%str(vaf*100)+'% Fail'

    def check_positive_control_VAF(self, row, cutoff=0.95):
        sample = row["Sample"]
        if not sample.endswith("SERACARE"):
            return "NA"

        total_rows = len(sera_care_df.index)
        var_rows = len(pos_control_variants.loc[(pos_control_variants['Sample'] == sample)
                            & pos_control_variants['Detected']].index)
        vaf = round(float(var_rows/total_rows),2)
        return "%s"%str(vaf*100)+'% Pass' if vaf >= cutoff else "%s"%str(vaf*100)+'% Fail'

    def build_report(self):
        doc_root = self.args.output
        global snp_overlap
        global unpaired_variants
        global hapmap_variants 
        global pos_control_variants
        global sera_care_df
        global excl_list
        excl_list = ["SERACARE", "HAPMAP", "NTC"]
        pattern = '|'.join(excl_list)
        # read data files
        mapping_reads = pd.read_csv(os.path.join(doc_root, "output/flagstat.tsv"), sep='\t')
        sample_mapping_reads = mapping_reads.copy()
        run_mapping_reads = mapping_reads.copy()
        sample_mapping_reads = sample_mapping_reads.loc[~(sample_mapping_reads['Sample'].str.contains(pattern, case=False))]
        run_mapping_reads = run_mapping_reads.loc[(run_mapping_reads['Sample'].str.contains(pattern, case=False))]

        dedup_reads = pd.read_csv(os.path.join(doc_root, "output/flagstat.dedup.tsv"), sep='\t')
        sample_dedup_reads = dedup_reads.copy()
        run_dedup_reads = dedup_reads.copy()
        sample_dedup_reads = sample_dedup_reads.loc[~(sample_dedup_reads['Sample'].str.contains(pattern, case=False))]
        run_dedup_reads = run_dedup_reads.loc[(run_dedup_reads['Sample'].str.contains(pattern, case=False))]

        coverage_df = pd.read_csv(os.path.join(doc_root, "output/coverage.samples.tsv"), sep='\t')
        sample_coverage_df = coverage_df.copy()
        run_coverage_df = coverage_df.copy()

        sample_coverage_df = sample_coverage_df.loc[~(sample_coverage_df['Sample'].str.contains(pattern, case=False))]
        sample_coverage_df = sample_coverage_df.loc[sample_coverage_df['Mode'] == "bed"]
        # print (sample_coverage_df.head(5).to_string())
        run_coverage_df = run_coverage_df.loc[(run_coverage_df['Sample'].str.contains(pattern, case=False))]
        run_coverage_df = run_coverage_df.loc[run_coverage_df['Mode'] == "bed"]

        snp_overlap = pd.read_csv(os.path.join(doc_root, "output/snp_overlap_80.tsv"), sep='\t')
        # filter on rows have overlap snps
        # snp_overlap = snp_overlap[snp_overlap.comb.str.contains('overlap')]

        # filter on variants for neg/pos controls
        unpaired_variants = pd.read_csv(os.path.join(doc_root, "output/annotations.unpaired.tsv"), sep='\t')
        pos_control_variants = unpaired_variants.copy()

        # filter on variants for neg/pos controls
        hapmap_variants = pd.read_csv(os.path.join(doc_root, "output/annotations/annotations.HapMap-Pool.tsv"), sep='\t') 

        # read tumor normal pair
        # sample_pairs = pd.read_csv(os.path.join(doc_root,"samples.tumor.normal.csv"),sep=',')
        sample_pairs = pd.read_csv(os.path.join(doc_root, "samples.tumor.normal.csv"), sep=',')
        sample_name_order = pd.concat([sample_pairs[['Normal']], sample_pairs[['Tumor']]], axis=1
                                      ).stack().reset_index(1, drop=True).to_frame('C').rename(index='CC{}'.format)
        print("reorder samples.....")
        print(sample_name_order.head(2).to_string())

        # SeraCare variants
        sera_care_df = pd.read_csv(self.args.SeraCare,sep='\t')
        sera_care_df['var_key'] = sera_care_df['CHROM'].str.cat(
            sera_care_df['POS'].astype(str), sep=":") \
            .str.cat(sera_care_df['REF'], sep="") \
            .str.cat(sera_care_df['ALT'], sep=">")

        unpaired_variants = unpaired_variants[(unpaired_variants['VariantCaller'] == "LoFreq")
                                              & (unpaired_variants['Func.refGene'] == "exonic")
                                              & (unpaired_variants['DP'] > 300)
                                              & (unpaired_variants['AF'] > 0.005)]

        hapmap_variants = hapmap_variants[(hapmap_variants['VariantCaller'] == "LoFreqSomatic")
                                              & (hapmap_variants['Func.refGene'] == "exonic") 
                                              & (hapmap_variants['DP'] > 200) 
                                              & (hapmap_variants['AF'] > 0.005)] 

        pos_control_variants = pos_control_variants[(pos_control_variants['VariantCaller'] == "LoFreq")
                                                    & (pos_control_variants['Func.refGene'] == "exonic")
                                                    & (pos_control_variants['Sample'].str.contains("SERACARE"))]
        pos_control_variants['var_key'] = pos_control_variants['CHROM'].str.cat(
            pos_control_variants['POS'].astype(str), sep=":") \
            .str.cat(pos_control_variants['REF'], sep="") \
            .str.cat(pos_control_variants['ALT'], sep=">")
        pos_control_variants = pos_control_variants.loc[pos_control_variants['var_key'].isin(sera_care_df['var_key'])]
        pos_control_variants['Detected'] = True
        print("pos controls: %s/%s" % (len(pos_control_variants.index), len(sera_care_df.index)))

        # check total mapping reads
        sample_mapping_reads = sample_mapping_reads[['Sample','TotalMappedReads']]
        run_mapping_reads = run_mapping_reads[['Sample','TotalMappedReads']]
        sample_mapping_reads['Mapped Reads QC'] = sample_mapping_reads.apply (lambda row: self.check_total_num_of_reads(row), axis=1)
        run_mapping_reads['Mapped Reads QC'] = run_mapping_reads.apply (lambda row: self.check_total_num_of_reads(row), axis=1)

        # check deduplicated reads
        sample_dedup_reads['Deduplicated Reads QC'] = sample_dedup_reads.apply (lambda row: self.check_total_num_of_reads(row), axis=1)
        run_dedup_reads['Deduplicated Reads QC'] = run_dedup_reads.apply (lambda row: self.check_total_num_of_reads(row), axis=1)

        sample_dedup_reads.rename(columns={"TotalMappedReads": "DeduplicatedReads"}, inplace=True)
        run_dedup_reads.rename(columns={"TotalMappedReads": "DeduplicatedReads"}, inplace=True)

        sample_dedup_reads = sample_dedup_reads[['Sample','DeduplicatedReads','Deduplicated Reads QC']]
        run_dedup_reads = run_dedup_reads[['Sample','DeduplicatedReads','Deduplicated Reads QC']]

        seq_QC_df = sample_mapping_reads.merge(sample_dedup_reads, on="Sample", how='inner')\
                    .merge(sample_coverage_df[['Sample','MeanCoverage','MedianCoverage','50']],
                                   on="Sample",how="inner")
        seq_QC_df['Deduplication_Rate'] = (1-seq_QC_df['DeduplicatedReads']/seq_QC_df['TotalMappedReads'])
        seq_QC_df['Deduplication Rate QC'] = seq_QC_df.apply (lambda row: self.check_dup_rate(row), axis=1)
        seq_QC_df = seq_QC_df[['Sample','TotalMappedReads','Mapped Reads QC','DeduplicatedReads','Deduplicated Reads QC','Deduplication Rate QC','MeanCoverage', 'MedianCoverage','50' ]]

        run_QC_df = run_mapping_reads.merge(run_dedup_reads, on="Sample", how='inner')\
                    .merge(run_coverage_df[['Sample','MeanCoverage','MedianCoverage','50']],
                                   on="Sample",how="inner")
        # check sample coverage
        seq_QC_df['Targets with <50 Coverage'] = seq_QC_df.apply (lambda row: self.check_sample_coverage(row), axis=1)
        # check mean and median coverage of sample
        seq_QC_df['Average and Median Coverage'] = seq_QC_df.apply (lambda row: self.check_average_median_coverage(row), axis=1)
        run_QC_df['Overall Coverage'] = run_QC_df.apply (lambda row: self.check_average_median_coverage(row), axis=1)
        # check overlapped homozygous SNPs
        seq_QC_df['Overlapped HOMO SNPs'] = seq_QC_df.apply (lambda row: self.check_snp_overlap(row), axis=1)
        # check private variants detected at coverage >200X and VAF >5%
        run_QC_df['Negative control QC'] = run_QC_df.apply (lambda row: self.check_negative_control_VAF(row), axis=1)
        # check >95% of these variants should have VAF within Mean± 2SD
        run_QC_df['Positive control QC'] = run_QC_df.apply (lambda row: self.check_positive_control_VAF(row), axis=1)
        # drop unused columns
        seq_QC_df.drop(columns=['TotalMappedReads','DeduplicatedReads','MedianCoverage','MeanCoverage','50'], inplace=True)
        seq_QC_df = seq_QC_df.drop_duplicates('Sample')
        run_QC_df.drop(columns=['TotalMappedReads','DeduplicatedReads','MedianCoverage','MeanCoverage','50'], inplace=True)
        run_QC_df = run_QC_df.drop_duplicates('Sample')
        # order sample by T->N
        seq_QC_df['name'] = seq_QC_df['Sample']
        seq_QC_df = seq_QC_df.set_index('name')
        seq_QC_df = seq_QC_df.loc[sample_name_order['C']]
        #seq_QC_df = seq_QC_df.loc[seq_QC_df['Sample'].isin(sample_name_order['C'])]
        print(seq_QC_df.head(2).to_string())
        seq_QC_df.insert(0, 'Row', (seq_QC_df.reset_index().index+1).astype(str))
        #seq_QC_df['Row'] = (seq_QC_df.reset_index().index+1).astype(str)
        pos_control_variants = pos_control_variants.merge(sera_care_df, on="var_key", how='right')
        pos_control_variants = pos_control_variants.fillna(0)
        pos_control_variants = pos_control_variants.round({"AF_x":2})
        pos_control_variants["AF 95% Interval (Historical Runs)"] = \
            round((pos_control_variants["AF.mean"]-2*pos_control_variants['AF.sd']),2).astype(str) + " ~ " + \
            round((pos_control_variants["AF.mean"]+2*pos_control_variants['AF.sd']),2).astype(str)
        pos_control_variants = pos_control_variants[['Gene','Coding', 'AAChange', 'DP', 'AF_x','AF 95% Interval (Historical Runs)','Detected']]
        pos_control_variants = pos_control_variants.rename({'AF_x':"AF"}, axis=1)
        #reset index
        pos_control_variants.insert(0, 'Row', (pos_control_variants.reset_index().index+1).astype(str))
        #rename sample
        run_QC_df["Sample"] = run_QC_df["Sample"].str.split("_", n=5, expand = True)[5]
        seq_QC_df["Sample"] = seq_QC_df["Sample"].str.split("_", n=5, expand = True)[5]
        seq_QC_df.to_csv(os.path.join(doc_root,'%s-QC.tsv'%self.args.project), sep="\t", index = False)

        return {"run_df":run_QC_df,"df":seq_QC_df, "SeraCare":pos_control_variants,"caseID":self.args.project,"run":self.args.run}
