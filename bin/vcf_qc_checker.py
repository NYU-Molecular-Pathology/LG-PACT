#!/usr/bin/env python3
"""
Complete QC Pipeline - Controls, Samples, HSMetrics, and Conpair QC Validation

Validates control samples first, then validates regular samples if controls pass.
Includes sample type classification, case-level QC, optional HSMetrics validation,
and optional Conpair concordance validation.
"""

import pandas as pd
import argparse
import sys
import re
from pathlib import Path
import csv


class ControlsQCChecker:
    """Validate controls against RunQC metrics."""
    
    def __init__(self, runqc_file, demux_file, output_file, runid=None):
        self.runqc_file = Path(runqc_file)
        self.demux_file = Path(demux_file)
        self.output_file = Path(output_file)
        self.runid = runid
        
        
    def log(self, message):
        """Logging disabled."""
        pass
    
    def parse_demux(self):
        """Parse Illumina demux sample sheet."""
        with open(self.demux_file, 'r') as f:
            lines = f.readlines()
        
        for i, line in enumerate(lines):
            if line.startswith('[Data]'):
                return pd.read_csv(self.demux_file, skiprows=i+1)
        
        raise ValueError("[Data] section not found in demux file")
    
    def is_control(self, sample_id):
        """Check if sample is a control."""
        s = str(sample_id)
        return ('_NTC_' in s or '_NC_' in s or '_SC_' in s or
                s.startswith('NTC_') or s.startswith('NC_') or s.startswith('SC_'))
    
    def extract_qc_status(self, qc_value, column_name=None):
        """Extract value and status from QC field."""
        val_str = str(qc_value).strip()
        
        if val_str in ['', 'nan', 'NA']:
            return 'NA', 'NA'
        elif 'Pass' in val_str:
            clean_val = val_str.replace(' Pass', '').strip()
            
            if column_name == 'Overall Coverage' and ',' in clean_val:
                clean_val = clean_val.split(',')[0].strip()
            
            return clean_val, 'Pass'
        elif 'Fail' in val_str:
            clean_val = val_str.replace(' Fail', '').strip()
            
            if column_name == 'Overall Coverage' and ',' in clean_val:
                clean_val = clean_val.split(',')[0].strip()
            
            return clean_val, 'Fail'
        else:
            if column_name == 'Overall Coverage' and ',' in val_str:
                val_str = val_str.split(',')[0].strip()
            
            return val_str, 'Unknown'
    
    def extract_control_info(self, sample_id):
        """Extract control name and type from Sample_ID."""
        parts = sample_id.split('_')
        
        for i, part in enumerate(parts):
            if part in ['NTC', 'NC', 'SC']:
                control_name = '_'.join(parts[i:])
                control_type = part
                return control_name, control_type
        
        return None, None
    
    def match_to_runqc(self, control_name, control_type, runqc_dict):
        """Match control to RunQC sample."""
        if control_name and control_name in runqc_dict:
            return runqc_dict[control_name], control_name
        
        for name in runqc_dict.keys():
            if control_type and control_type in name:
                return runqc_dict[name], name
        
        return None, None
    
    def process_control(self, ctrl_row, runqc_dict, runqc_columns):
        """Process single control sample."""
        sample_id = str(ctrl_row['Sample_ID'])
        case_id = str(ctrl_row['Case_ID'])
        
        control_name, control_type = self.extract_control_info(sample_id)
        runqc_match, matched_name = self.match_to_runqc(control_name, control_type, runqc_dict)
        
        result = {
            'Sample_ID': sample_id,
            'Case_ID': case_id,
            'Control_Type': control_type or 'Unknown',
            'RunQC_Sample': matched_name or 'Not Found',
            'Matched': 'Yes' if runqc_match else 'No'
        }
        
        # Add Run_ID if provided
        if self.runid:
            result['Run_ID'] = self.runid
        
        if runqc_match:
            failed_metrics = []
            
            for col in runqc_columns:
                if col != 'Sample':
                    val, status = self.extract_qc_status(runqc_match[col], column_name=col)
                    result[f'{col}_Value'] = val
                    result[f'{col}_Status'] = status
                    
                    if status == 'Fail':
                        failed_metrics.append(col)
            
            all_statuses = [v for k, v in result.items() if k.endswith('_Status')]
            
            if 'Fail' in all_statuses:
                result['Overall_QC'] = 'FAIL'
                result['Failed_Metrics'] = '; '.join(failed_metrics)
            elif all(s == 'NA' for s in all_statuses):
                result['Overall_QC'] = 'NO_DATA'
                result['Failed_Metrics'] = ''
            else:
                result['Overall_QC'] = 'PASS'
                result['Failed_Metrics'] = ''
        else:
            result['Overall_QC'] = 'NOT_FOUND'
            result['Failed_Metrics'] = ''
        
        return result
    
    def run(self):
        """Execute QC check and return results."""
        self.log(f"Reading demux file: {self.demux_file}")
        demux_df = self.parse_demux()
        
        self.log(f"Reading RunQC file: {self.runqc_file}")
        runqc_df = pd.read_csv(self.runqc_file, sep='\t')
        
        control_mask = demux_df['Sample_ID'].apply(self.is_control)
        controls_df = demux_df[control_mask]
        
        if controls_df.empty:
            print("ERROR: No controls found in demux file", file=sys.stderr)
            return None
        
        self.log(f"Found {len(controls_df)} controls")
        
        runqc_dict = {row['Sample']: row.to_dict() for _, row in runqc_df.iterrows()}
        
        results = []
        for _, ctrl_row in controls_df.iterrows():
            result = self.process_control(ctrl_row, runqc_dict, runqc_df.columns)
            results.append(result)
        
        results_df = pd.DataFrame(results)
        
        self.output_file.parent.mkdir(parents=True, exist_ok=True)
        results_df.to_csv(self.output_file, index=False)
        
        print(f"Controls QC saved: {self.output_file}")
        self.print_summary(results_df)
        
        return results_df
    
    def print_summary(self, results_df):
        """Print QC summary."""
        print("\n" + "="*70)
        print("CONTROLS QC SUMMARY")
        print("="*70)
        
        for _, row in results_df.iterrows():
            status = row['Overall_QC']
            print(f"{row['Control_Type']:3s} | {row['RunQC_Sample']:15s} | {status}")
            
            if row['Failed_Metrics']:
                print(f"    Failed: {row['Failed_Metrics']}")
        
        print("="*70)


class SamplesQCChecker:
    """Validate regular samples against QC metrics."""
    
    def __init__(self, qc_file, demux_file, output_file, runid=None):
        self.qc_file = Path(qc_file)
        self.demux_file = Path(demux_file)
        self.output_file = Path(output_file)
        self.runid = runid
        
        
    def log(self, message):
        """Logging disabled."""
        pass
    
    def parse_demux(self):
        """Parse Illumina demux sample sheet."""
        with open(self.demux_file, 'r') as f:
            lines = f.readlines()
        
        for i, line in enumerate(lines):
            if line.startswith('[Data]'):
                return pd.read_csv(self.demux_file, skiprows=i+1)
        
        raise ValueError("[Data] section not found in demux file")
    
    def is_sample(self, sample_id):
        """Check if sample is a regular sample (not control)."""
        s = str(sample_id)
        is_control = ('_NTC_' in s or '_NC_' in s or '_SC_' in s or
                     s.startswith('NTC_') or s.startswith('NC_') or s.startswith('SC_'))
        return not is_control
    
    def parse_coverage_column(self, coverage_str):
        """Parse coverage column with comma-separated average,median values.
        
        Expected format: "average,median Pass/Fail" or "average,median"
        Example: "790.55,603 Pass" or "790.55,603"
        Returns: (average, median, status)
        """
        val_str = str(coverage_str).strip()
        
        if val_str in ['', 'nan', 'NA']:
            return 'NA', 'NA', 'NA'
        
        # Extract status
        status = 'Unknown'
        if 'Pass' in val_str:
            status = 'Pass'
            val_str = val_str.replace(' Pass', '').strip()
        elif 'Fail' in val_str:
            status = 'Fail'
            val_str = val_str.replace(' Fail', '').strip()
        
        # Split on comma to get average and median
        if ',' in val_str:
            parts = val_str.split(',')
            average = parts[0].strip()
            median = parts[1].strip() if len(parts) > 1 else 'NA'
            return average, median, status
        else:
            # Only one value provided (treat as average)
            return val_str, 'NA', status
    
    def extract_qc_status(self, qc_value, column_name=None):
        """Extract value and status from QC field."""
        val_str = str(qc_value).strip()
        
        if val_str in ['', 'nan', 'NA']:
            return 'NA', 'NA'
        elif 'Pass' in val_str:
            clean_val = val_str.replace(' Pass', '').strip()
            
            if column_name == 'Overall Coverage' and ',' in clean_val:
                clean_val = clean_val.split(',')[0].strip()
            
            return clean_val, 'Pass'
        elif 'Fail' in val_str:
            clean_val = val_str.replace(' Fail', '').strip()
            
            if column_name == 'Overall Coverage' and ',' in clean_val:
                clean_val = clean_val.split(',')[0].strip()
            
            return clean_val, 'Fail'
        else:
            if column_name == 'Overall Coverage' and ',' in val_str:
                val_str = val_str.split(',')[0].strip()
            
            return val_str, 'Unknown'
    
    def extract_case_id(self, sample_id):
        """Extract case ID from Sample_ID."""
        parts = sample_id.split('_')
        
        for i, part in enumerate(parts):
            if part.startswith(('TS', 'WS', 'BS')):
                return '_'.join(parts[i:])
        
        return None
    
    def match_to_qc(self, case_id, qc_dict):
        """Match sample to QC file."""
        if case_id in qc_dict:
            return qc_dict[case_id], case_id
        
        for qc_sample_id in qc_dict.keys():
            if case_id in qc_sample_id or qc_sample_id in case_id:
                return qc_dict[qc_sample_id], qc_sample_id
        
        return None, None
    
    def classify_sample_type(self, sample_row, qc_match):
        """Classify sample as Tumor or Normal based on Case_ID and Tumor_Content."""
        case_id = str(sample_row['Case_ID'])
        
        if 'Tumor_Content' in sample_row:
            tumor_content = str(sample_row['Tumor_Content'])
            
            if tumor_content != '0' and tumor_content != '0.0' and tumor_content != 'nan':
                try:
                    tc_value = float(tumor_content)
                    if tc_value > 0:
                        return 'Tumor'
                except (ValueError, TypeError):
                    pass
        
        if '-A' in case_id or '-B' in case_id or '-C' in case_id or '-D' in case_id:
            return 'Tumor'
        
        return 'Normal'
    
    def perform_case_level_qc(self, samples_df, qc_dict):
        """Perform case-level QC: Check that each case has both tumor and normal."""
        if 'TM_Number' not in samples_df.columns:
            return {}
        
        case_qc = {}
        
        for test_num in samples_df['TM_Number'].unique():
            if pd.isna(test_num):
                continue
            
            test_samples = samples_df[samples_df['TM_Number'] == test_num]
            
            sample_types = []
            for _, sample_row in test_samples.iterrows():
                case_id = str(sample_row['Case_ID'])
                qc_match, _ = self.match_to_qc(case_id, qc_dict)
                sample_type = self.classify_sample_type(sample_row, qc_match)
                sample_types.append(sample_type)
            
            has_tumor = 'Tumor' in sample_types
            has_normal = 'Normal' in sample_types
            
            if has_tumor and has_normal:
                case_qc[test_num] = 'PASS'
            else:
                case_qc[test_num] = 'FAIL'
        
        return case_qc
    
    def process_sample(self, sample_row, qc_dict, qc_columns, case_qc_dict):
        """Process single sample."""
        sample_id = str(sample_row['Sample_ID'])
        case_id = str(sample_row['Case_ID'])
        case_part = self.extract_case_id(sample_id)
        qc_match, matched_name = self.match_to_qc(case_id, qc_dict)
        
        sample_type = self.classify_sample_type(sample_row, qc_match)
        
        result = {
            'Sample_ID': sample_id,
            'Case_ID': case_id,
            'Sample_Type': sample_type,
            'QC_Sample': matched_name or 'Not Found',
            'Matched': 'Yes' if qc_match else 'No'
        }
        
        # Add Run_ID if provided
        if self.runid:
            result['Run_ID'] = self.runid
        
        if 'EPIC_ID' in sample_row:
            result['EPIC_ID'] = str(sample_row['EPIC_ID'])
        if 'TM_Number' in sample_row:
            test_num = str(sample_row['TM_Number'])
            result['TM_Number'] = test_num
        
        if 'Tumor_Content' in sample_row:
            result['Tumor_Content'] = str(sample_row['Tumor_Content'])
        
        if qc_match:
            failed_metrics = []
            
            for col in qc_columns:
                if col != 'Sample':
                    # Special handling for coverage columns with comma-separated values
                    if col in ['Average Coverage', 'Median Coverage', 'Average and Median Coverage']:
                        average, median, status = self.parse_coverage_column(qc_match[col])
                        
                        # For "Average and Median Coverage", create clearer column names
                        if col == 'Average and Median Coverage':
                            result['Average_Coverage'] = average
                            result['Median_Coverage'] = median
                            result[f'{col}_Status'] = status
                        else:
                            result[f'{col}_Average_Value'] = average
                            result[f'{col}_Median_Value'] = median
                            result[f'{col}_Status'] = status
                        
                        if status == 'Fail':
                            failed_metrics.append(col)
                    else:
                        val, status = self.extract_qc_status(qc_match[col], column_name=col)
                        result[f'{col}_Value'] = val
                        result[f'{col}_Status'] = status
                        
                        if status == 'Fail':
                            failed_metrics.append(col)
            
            all_statuses = [v for k, v in result.items() if k.endswith('_Status') and k != 'Case_Level_QC']
            
            if 'Fail' in all_statuses:
                result['Sample_Level_QC'] = 'FAIL'
                result['Failed_Metrics'] = '; '.join(failed_metrics)
            elif all(s == 'NA' for s in all_statuses):
                result['Sample_Level_QC'] = 'NO_DATA'
                result['Failed_Metrics'] = ''
            else:
                result['Sample_Level_QC'] = 'PASS'
                result['Failed_Metrics'] = ''
        else:
            result['Sample_Level_QC'] = 'NOT_FOUND'
            result['Failed_Metrics'] = ''
        
        return result
    
    def run(self):
        """Execute QC check and return results."""
        self.log(f"Reading demux file: {self.demux_file}")
        demux_df = self.parse_demux()
        
        self.log(f"Reading QC file: {self.qc_file}")
        qc_df = pd.read_csv(self.qc_file, sep='\t')
        
        sample_mask = demux_df['Sample_ID'].apply(self.is_sample)
        samples_df = demux_df[sample_mask]
        
        if samples_df.empty:
            print("WARNING: No regular samples found in demux file")
            return None
        
        self.log(f"Found {len(samples_df)} samples")
        
        qc_dict = {row['Sample']: row.to_dict() for _, row in qc_df.iterrows()}
        
        case_qc_dict = self.perform_case_level_qc(samples_df, qc_dict)
        
        results = []
        for _, sample_row in samples_df.iterrows():
            result = self.process_sample(sample_row, qc_dict, qc_df.columns, case_qc_dict)
            results.append(result)
        
        results_df = pd.DataFrame(results)
        
        self.output_file.parent.mkdir(parents=True, exist_ok=True)
        results_df.to_csv(self.output_file, index=False)
        
        print(f"Samples QC saved: {self.output_file}")
        self.print_summary(results_df)
        
        return results_df
    
    def print_summary(self, results_df):
        """Print QC summary."""
        print("\n" + "="*70)
        print("SAMPLES QC SUMMARY")
        print("="*70)
        
        for _, row in results_df.iterrows():
            specimen = row['Case_ID']
            sample_type = row['Sample_Type']
            sample_qc = row.get('Sample_Level_QC', 'NA')
            
            print(f"{specimen:20s} | {sample_type:6s} | Sample QC: {sample_qc}")
            
            if row.get('Failed_Metrics', ''):
                print(f"    Failed: {row['Failed_Metrics']}")
        
        print("="*70)


class HSMetricsQCChecker:
    """Validate samples against HSMetrics coverage metrics."""
    
    def __init__(self, hsmetrics_file, demux_file, output_file, threshold=0.75, runid=None):
        self.hsmetrics_file = Path(hsmetrics_file)
        self.demux_file = Path(demux_file)
        self.output_file = Path(output_file)
        self.threshold = threshold
        self.runid = runid
        
        
    def log(self, message):
        """Logging disabled."""
        pass
    
    def parse_demux(self):
        """Parse Illumina demux sample sheet."""
        with open(self.demux_file, 'r') as f:
            lines = f.readlines()
        
        for i, line in enumerate(lines):
            if line.startswith('[Data]'):
                return pd.read_csv(self.demux_file, skiprows=i+1)
        
        raise ValueError("[Data] section not found in demux file")
    
    def is_sample(self, sample_id):
        """Check if sample is a regular sample (not control)."""
        s = str(sample_id)
        is_control = ('_NTC_' in s or '_NC_' in s or '_SC_' in s or
                     s.startswith('NTC_') or s.startswith('NC_') or s.startswith('SC_'))
        return not is_control
    
    def extract_case_id(self, sample_id):
        """Extract case ID from Sample_ID."""
        parts = sample_id.split('_')
        
        for i, part in enumerate(parts):
            if part.startswith(('TS', 'WS', 'BS')):
                return '_'.join(parts[i:])
        
        return None
    
    def match_to_hsmetrics(self, case_id, hsmetrics_dict):
        """Match sample to HSMetrics file."""
        if case_id in hsmetrics_dict:
            return hsmetrics_dict[case_id], case_id
        
        for hs_sample_id in hsmetrics_dict.keys():
            if case_id in hs_sample_id or hs_sample_id in case_id:
                return hsmetrics_dict[hs_sample_id], hs_sample_id
        
        return None, None
    
    def process_sample(self, sample_row, hsmetrics_dict):
        """Process single sample."""
        sample_id = str(sample_row['Sample_ID'])
        case_id = str(sample_row['Case_ID'])
        case_part = self.extract_case_id(sample_id)
        hs_match, matched_name = self.match_to_hsmetrics(case_id, hsmetrics_dict)
        
        result = {
            'Sample_ID': sample_id,
            'Case_ID': case_id,
            'HSMetrics_Sample': matched_name or 'Not Found',
            'Matched': 'Yes' if hs_match else 'No'
        }
        
        # Add Run_ID if provided
        if self.runid:
            result['Run_ID'] = self.runid
        
        if 'EPIC_ID' in sample_row:
            result['EPIC_ID'] = str(sample_row['EPIC_ID'])
        if 'TM_Number' in sample_row:
            result['TM_Number'] = str(sample_row['TM_Number'])
        if 'Tumor_Content' in sample_row:
            result['Tumor_Content'] = str(sample_row['Tumor_Content'])
        
        if hs_match:
            pct_250x = hs_match.get('PCT_TARGET_BASES_250X', None)
            
            if pct_250x is not None and pct_250x != '':
                try:
                    pct_value = float(pct_250x)
                    result['PCT_TARGET_BASES_250X'] = pct_value
                    
                    if pct_value >= self.threshold:
                        result['QC_Status'] = 'PASS'
                    else:
                        result['QC_Status'] = 'FAIL'
                    
                except (ValueError, TypeError):
                    result['PCT_TARGET_BASES_250X'] = pct_250x
                    result['QC_Status'] = 'UNKNOWN'
            else:
                result['PCT_TARGET_BASES_250X'] = 'NA'
                result['QC_Status'] = 'NO_DATA'
        else:
            result['PCT_TARGET_BASES_250X'] = 'NA'
            result['QC_Status'] = 'NOT_FOUND'
        
        return result
    
    def run(self):
        """Execute HSMetrics QC check."""
        self.log(f"Reading demux file: {self.demux_file}")
        demux_df = self.parse_demux()
        
        self.log(f"Reading HSMetrics file: {self.hsmetrics_file}")
        hsmetrics_df = pd.read_csv(self.hsmetrics_file)
        
        sample_mask = demux_df['Sample_ID'].apply(self.is_sample)
        samples_df = demux_df[sample_mask]
        
        if samples_df.empty:
            print("WARNING: No regular samples found in demux file")
            return None
        
        self.log(f"Found {len(samples_df)} samples")
        
        hsmetrics_dict = {}
        for _, row in hsmetrics_df.iterrows():
            sample_id = row.get('SampleID', None)
            if sample_id:
                hsmetrics_dict[sample_id] = row.to_dict()
        
        self.log(f"Found {len(hsmetrics_dict)} samples in HSMetrics file")
        
        results = []
        for _, sample_row in samples_df.iterrows():
            result = self.process_sample(sample_row, hsmetrics_dict)
            results.append(result)
        
        results_df = pd.DataFrame(results)
        
        self.output_file.parent.mkdir(parents=True, exist_ok=True)
        results_df.to_csv(self.output_file, index=False)
        
        print(f"HSMetrics QC saved: {self.output_file}")
        self.print_summary(results_df)
        
        return results_df
    
    def print_summary(self, results_df):
        """Print QC summary."""
        print("\n" + "="*70)
        print("HSMETRICS QC SUMMARY")
        print(f"Threshold: PCT_TARGET_BASES_250X >= {self.threshold:.0%}")
        print("="*70)
        
        for _, row in results_df.iterrows():
            specimen = row['Case_ID']
            pct_250x = row['PCT_TARGET_BASES_250X']
            status = row['QC_Status']
            
            if isinstance(pct_250x, (int, float)):
                print(f"{specimen:20s} | {pct_250x:.4f} ({pct_250x*100:.2f}%) | {status}")
            else:
                print(f"{specimen:20s} | {pct_250x:>20s} | {status}")
        
        print("="*70)
        
        pass_count = len(results_df[results_df['QC_Status'] == 'PASS'])
        fail_count = len(results_df[results_df['QC_Status'] == 'FAIL'])
        total = len(results_df)
        
        print(f"\nTotal: {total} | Pass: {pass_count} | Fail: {fail_count}")


class ConpairQCChecker:
    """Validate samples against Conpair concordance metrics."""
    
    def __init__(self, conpair_file, demux_file, output_file, 
                 concordance_threshold=75.0, runid=None):
        self.conpair_file = Path(conpair_file)
        self.demux_file = Path(demux_file)
        self.output_file = Path(output_file)
        self.concordance_threshold = concordance_threshold
        self.runid = runid
        
        
    def log(self, message):
        """Logging disabled."""
        pass
    
    def parse_demux(self):
        """Parse Illumina demux sample sheet."""
        with open(self.demux_file, 'r') as f:
            lines = f.readlines()
        
        for i, line in enumerate(lines):
            if line.startswith('[Data]'):
                return pd.read_csv(self.demux_file, skiprows=i+1)
        
        raise ValueError("[Data] section not found in demux file")
    
    def is_sample(self, sample_id):
        """Check if sample is a regular sample (not control)."""
        s = str(sample_id)
        is_control = ('_NTC_' in s or '_NC_' in s or '_SC_' in s or
                     s.startswith('NTC_') or s.startswith('NC_') or s.startswith('SC_'))
        return not is_control
    
    def is_tumor_sample(self, case_id):
        """Check if sample is a tumor sample based on Case_ID."""
        spec_str = str(case_id)
        # Tumor samples have suffix like -A2, -B1, -C2, -D2 (letter + number)
        return bool(re.search(r'-[ABCD]\d+', spec_str))
    
    def match_to_conpair(self, sample_id, case_id, conpair_dict):
        """Match sample to Conpair file using multiple strategies."""
        for conpair_id, conpair_data in conpair_dict.items():
            # Strategy 1: Check if conpair_id is in sample_id (handles suffixes)
            if conpair_id in sample_id:
                return conpair_data, conpair_id
            
            # Strategy 2: Check if sample_id is in conpair_id (handles concatenated IDs)
            if sample_id in conpair_id:
                return conpair_data, conpair_id
            
            # Strategy 3: Check if case_id is in conpair_id (handles fragmented IDs)
            if case_id in conpair_id:
                return conpair_data, conpair_id
        
        return None, None
    
    def process_sample(self, sample_row, conpair_dict):
        """Process a single sample."""
        sample_id = str(sample_row['Sample_ID'])
        case_id = str(sample_row['Case_ID'])
        
        is_tumor = self.is_tumor_sample(case_id)
        sample_type = 'Tumor' if is_tumor else 'Normal'
        
        result = {
            'Sample_ID': sample_id,
            'Case_ID': case_id,
            'Sample_Type': sample_type,
            'Conpair_Sample': 'Not Found',
            'Matched': 'No'
        }
        
        # Add Run_ID if provided
        if self.runid:
            result['Run_ID'] = self.runid
        
        # Add optional columns
        if 'EPIC_ID' in sample_row:
            result['EPIC_ID'] = str(sample_row['EPIC_ID'])
        if 'TM_Number' in sample_row:
            result['TM_Number'] = str(sample_row['TM_Number'])
        if 'Tumor_Content' in sample_row:
            result['Tumor_Content'] = str(sample_row['Tumor_Content'])
        
        # Match to conpair using Sample_ID and Case_ID
        conpair_match, matched_id = self.match_to_conpair(sample_id, case_id, conpair_dict)
        
        if conpair_match:
            result['Conpair_Sample'] = matched_id
            result['Matched'] = 'Yes'
            
            # Parse concordance
            concordance_str = str(conpair_match.get('Conpair_Concordance', '')).strip()
            try:
                concordance_val = float(concordance_str)
            except (ValueError, TypeError):
                concordance_val = None
            
            # Store concordance
            result['Concordance_Percent'] = concordance_val if concordance_val is not None else 'NA'
            
            # Concordance check
            if concordance_val is not None:
                if concordance_val > self.concordance_threshold:
                    result['Concordance_Status'] = 'PASS'
                else:
                    result['Concordance_Status'] = 'FAIL'
            else:
                result['Concordance_Status'] = 'NO_DATA'
            
            # Overall QC is same as Concordance Status
            result['Overall_QC'] = result['Concordance_Status']
            
            # Failed metrics
            failed = []
            if result['Concordance_Status'] == 'FAIL':
                failed.append('Concordance')
            
            result['Failed_Metrics'] = '; '.join(failed) if failed else ''
            
        else:
            result['Concordance_Percent'] = 'NA'
            result['Concordance_Status'] = 'NOT_FOUND'
            result['Overall_QC'] = 'NOT_FOUND'
            result['Failed_Metrics'] = ''
        
        return result
    
    def run(self):
        """Execute Conpair QC check."""
        self.log(f"Reading demux file: {self.demux_file}")
        demux_df = self.parse_demux()
        
        self.log(f"Reading Conpair file: {self.conpair_file}")
        conpair_df = pd.read_csv(self.conpair_file)
        
        # Filter for samples only (not controls)
        sample_mask = demux_df['Sample_ID'].apply(self.is_sample)
        samples_df = demux_df[sample_mask]
        
        if samples_df.empty:
            print("WARNING: No regular samples found in demux file")
            return None
        
        self.log(f"Found {len(samples_df)} samples")
        
        # Create Conpair lookup
        conpair_dict = {}
        for _, row in conpair_df.iterrows():
            sample_id = row.get('SampleID', None)
            if sample_id:
                conpair_dict[sample_id] = row.to_dict()
        
        self.log(f"Found {len(conpair_dict)} entries in Conpair file")
        
        # Process each sample
        results = []
        for _, sample_row in samples_df.iterrows():
            result = self.process_sample(sample_row, conpair_dict)
            results.append(result)
        
        if not results:
            print("WARNING: No samples could be processed")
            return None
        
        results_df = pd.DataFrame(results)
        
        # Save output
        self.output_file.parent.mkdir(parents=True, exist_ok=True)
        results_df.to_csv(self.output_file, index=False)
        
        print(f"Conpair QC saved: {self.output_file}")
        self.print_summary(results_df)
        
        return results_df
    
    def print_summary(self, results_df):
        """Print QC summary."""
        print("\n" + "="*70)
        print("CONPAIR QC SUMMARY")
        print(f"Concordance Threshold: > {self.concordance_threshold}%")
        print("="*70)
        
        for _, row in results_df.iterrows():
            specimen = row['Case_ID']
            sample_type = row['Sample_Type']
            concordance = row['Concordance_Percent']
            status = row['Overall_QC']
            
            if isinstance(concordance, (int, float)):
                print(f"{specimen:20s} | {sample_type:6s} | Concordance: {concordance:.2f}% | {status}")
            else:
                print(f"{specimen:20s} | {sample_type:6s} | Concordance: {concordance:>10s} | {status}")
            
            if row['Failed_Metrics']:
                print(f"    Failed: {row['Failed_Metrics']}")
        
        print("="*70)
        
        # Overall stats
        pass_count = len(results_df[results_df['Overall_QC'] == 'PASS'])
        fail_count = len(results_df[results_df['Overall_QC'] == 'FAIL'])
        not_found_count = len(results_df[results_df['Overall_QC'] == 'NOT_FOUND'])
        total = len(results_df)
        
        print(f"\nTotal: {total} | Pass: {pass_count} | Fail: {fail_count} | Not Found: {not_found_count}")

class ExomeCoverageQCChecker:
    """Check exome coverage percentages against threshold (HSMetrics-style format)."""
    
    def __init__(self, tumor_file, normal_file, demux_file, output_file, threshold=90.0, runid=None):
        self.tumor_file = Path(tumor_file)
        self.normal_file = Path(normal_file)
        self.demux_file = Path(demux_file)
        self.output_file = Path(output_file)
        self.threshold = threshold
        self.runid = runid
        
    
    def log(self, message):
        """Logging disabled."""
        pass
    
    def is_sample(self, sample_id):
        """Check if sample is a regular sample (not control)."""
        s = str(sample_id)
        is_control = ('_NTC_' in s or '_NC_' in s or '_SC_' in s or
                     s.startswith('NTC_') or s.startswith('NC_') or s.startswith('SC_'))
        return not is_control
    
    def classify_sample_type(self, sample_row):
        """Classify sample as Tumor or Normal based on Case_ID and Tumor_Content."""
        case_id = str(sample_row.get('Case_ID', ''))
        
        # Check Tumor_Content if available
        if 'Tumor_Content' in sample_row:
            tumor_content = str(sample_row['Tumor_Content'])
            
            if tumor_content not in ('0', '0.0', 'nan', '', 'NA'):
                try:
                    tc_value = float(tumor_content)
                    if tc_value > 0:
                        return 'Tumor'
                except (ValueError, TypeError):
                    pass
        
        # Check for tumor suffixes (-A, -B, -C, -D)
        if '-A' in case_id or '-B' in case_id or '-C' in case_id or '-D' in case_id:
            return 'Tumor'
        
        return 'Normal'
    
    def parse_demux(self):
        """Parse demux sample sheet using readlines - similar to HSMetricsQCChecker."""
        with open(self.demux_file, 'r') as f:
            lines = f.readlines()
        
        # Look for [Data] section (Illumina format)
        data_start = None
        for i, line in enumerate(lines):
            if line.strip().startswith('[Data]'):
                data_start = i + 1
                break
        
        # If no [Data] section, assume it's a simple CSV starting at line 0
        if data_start is None:
            data_start = 0
        
        # Parse CSV data using csv.DictReader
        demux_samples = []
        csv_reader = csv.DictReader(lines[data_start:])
        
        for row in csv_reader:
            demux_samples.append(row)
        
        self.log(f"Parsed {len(demux_samples)} samples from demux")
        if False and len(demux_samples) > 0:
            print(f"  Demux columns: {list(demux_samples[0].keys())}")
            print(f"  Sample demux Sample_IDs:")
            for sample in demux_samples[:3]:
                print(f"    - {sample.get('Sample_ID', 'N/A')}")
        
        return demux_samples
    
    def read_coverage_file(self, coverage_file, sample_type):
        """Read coverage file using readlines and csv module."""
        if not coverage_file.exists():
            self.log(f"WARNING: Coverage file not found: {coverage_file}")
            return {}
        
        coverage_dict = {}
        
        try:
            with open(coverage_file, 'r') as f:
                lines = f.readlines()
            
            # Parse CSV using csv.DictReader
            csv_reader = csv.DictReader(lines)
            
            # Show columns found
            fieldnames = csv_reader.fieldnames
            self.log(f"{sample_type} file columns: {fieldnames}")
            
            # Determine which columns to use
            sample_id_col = None
            coverage_col = None
            
            if 'SampleID' in fieldnames:
                sample_id_col = 'SampleID'
            elif 'Sample_ID' in fieldnames:
                sample_id_col = 'Sample_ID'
            
            if 'PercentCoverage' in fieldnames:
                coverage_col = 'PercentCoverage'
            elif 'Exome_Percent_Coverage' in fieldnames:
                coverage_col = 'Exome_Percent_Coverage'
            
            if not sample_id_col:
                print(f"  ❌ ERROR: No Sample ID column found in {sample_type} file!")
                print(f"  Available columns: {fieldnames}")
                print(f"  Expected: 'SampleID' or 'Sample_ID'")
                return {}
            
            if not coverage_col:
                print(f"  ❌ ERROR: No coverage column found in {sample_type} file!")
                print(f"  Available columns: {fieldnames}")
                print(f"  Expected: 'PercentCoverage' or 'Exome_Percent_Coverage'")
                return {}
            
            self.log(f"Using columns: SampleID='{sample_id_col}', Coverage='{coverage_col}'")
            
            count = 0
            for row in csv_reader:
                # Get sample ID using detected column
                sample_id = row.get(sample_id_col, '').strip()
                
                # Get coverage percentage using detected column
                pct_cov_str = row.get(coverage_col, '').strip()
                
                # Debug: show first few rows
                if False and count < 3:
                    print(f"    Row {count+1}: {sample_id_col}='{sample_id}', {coverage_col}='{pct_cov_str}'")
                
                if sample_id:
                    # Store the value (keep as string, will convert later)
                    pct_cov = pct_cov_str if pct_cov_str else None
                    coverage_dict[sample_id] = {
                        'Exome_Percent_Coverage': pct_cov,
                        'Sample_Type': sample_type
                    }
                    count += 1
            
            self.log(f"Read {count} {sample_type} samples from {coverage_file.name}")
            
            # Show sample IDs if verbose
            if False and count > 0:
                print(f"  Sample {sample_type.lower()} coverage IDs:")
                for i, sid in enumerate(list(coverage_dict.keys())[:3]):
                    cov_val = coverage_dict[sid]['Exome_Percent_Coverage']
                    print(f"    - {sid} → Coverage: {cov_val}")
            
        except Exception as e:
            print(f"WARNING: Failed to read {sample_type} coverage file: {e}")
            import traceback
            traceback.print_exc()
        
        return coverage_dict
    
    def read_coverage_files(self):
        """Read and return tumor and normal coverage files as separate dictionaries."""
        # Read tumor file
        self.log(f"Reading tumor coverage file...")
        tumor_dict = self.read_coverage_file(self.tumor_file, 'Tumor')
        
        # Read normal file
        self.log(f"Reading normal coverage file...")
        normal_dict = self.read_coverage_file(self.normal_file, 'Normal')
        
        total_samples = len(tumor_dict) + len(normal_dict)
        self.log(f"Total coverage samples loaded: {total_samples} (Tumor: {len(tumor_dict)}, Normal: {len(normal_dict)})")
        
        return tumor_dict, normal_dict
    
    def match_to_coverage(self, sample_id, case_id, sample_type, tumor_dict, normal_dict):
        """
        Match demux sample to coverage file using string contains.
        Matches tumor samples to tumor coverage, normal samples to normal coverage.
        
        Primary matching: Check if Case_ID from demux is contained in 
        SampleID from exome coverage data.
        
        Example:
          Demux Case_ID: "TS25-12345-A2" (Tumor)
          Tumor Coverage SampleID: "TS25-12345-A2_B25-1122"
          Check: "TS25-12345-A2" in "TS25-12345-A2_B25-1122" → Match!
        """
        
        # Select the appropriate coverage dictionary based on sample type
        if sample_type == 'Tumor':
            coverage_dict = tumor_dict
            self.log(f"Matching {case_id} to TUMOR coverage")
        else:
            coverage_dict = normal_dict
            self.log(f"Matching {case_id} to NORMAL coverage")
        
        # Step 1: Try exact match with Case_ID
        if case_id in coverage_dict:
            self.log(f"Exact match: Case_ID '{case_id}' found in {sample_type} coverage")
            return coverage_dict[case_id], case_id
        
        # Step 2: Try exact match with Sample_ID
        if sample_id in coverage_dict:
            self.log(f"Exact match: Sample_ID '{sample_id}' found in {sample_type} coverage")
            return coverage_dict[sample_id], sample_id
        
        # Step 3: PRIMARY - String contains matching with Case_ID
        # Check if demux Case_ID is contained in coverage SampleID
        for cov_sample_id, cov_data in coverage_dict.items():
            if case_id in cov_sample_id:
                self.log(f"String contains match: Case_ID '{case_id}' IN {sample_type} coverage '{cov_sample_id}'")
                return cov_data, cov_sample_id
        
        # Step 4: String contains matching with Sample_ID
        # Check if demux Sample_ID is contained in coverage SampleID
        for cov_sample_id, cov_data in coverage_dict.items():
            if sample_id in cov_sample_id:
                self.log(f"String contains match: Sample_ID '{sample_id}' IN {sample_type} coverage '{cov_sample_id}'")
                return cov_data, cov_sample_id
        
        # Step 5: Reverse contains (fallback)
        # Check if coverage SampleID is contained in demux IDs
        for cov_sample_id, cov_data in coverage_dict.items():
            if cov_sample_id in case_id or cov_sample_id in sample_id:
                self.log(f"Reverse contains match: {sample_type} coverage '{cov_sample_id}' IN demux IDs")
                return cov_data, cov_sample_id
        
        # No match found
        self.log(f"No match found in {sample_type} coverage for Sample_ID: '{sample_id}', Case_ID: '{case_id}'")
        return None, None
    
    def process_sample(self, sample_row, tumor_dict, normal_dict):
        """Process single sample from demux file."""
        sample_id = sample_row.get('Sample_ID', '')
        case_id = sample_row.get('Case_ID', sample_id)
        
        # Classify sample type (Tumor or Normal)
        sample_type = self.classify_sample_type(sample_row)
        self.log(f"Classified {case_id} as {sample_type}")
        
        # Match to coverage data based on sample type
        cov_match, matched_name = self.match_to_coverage(sample_id, case_id, sample_type, tumor_dict, normal_dict)
        
        result = {
            'Sample_ID': sample_id,
            'Case_ID': case_id,
            'Sample_Type': sample_type,  # Add detected sample type
            'Exome_Sample': matched_name if matched_name else 'Not Found',
            'Matched': 'Yes' if cov_match else 'No'
        }
        
        # Add Run_ID if provided
        if self.runid:
            result['Run_ID'] = self.runid
        
        # Add optional demux columns
        if 'EPIC_ID' in sample_row:
            result['EPIC_ID'] = sample_row['EPIC_ID']
        if 'TM_Number' in sample_row:
            result['TM_Number'] = sample_row['TM_Number']
        # Don't override Sample_Type if already set from demux
        # We prefer our detected type, but keep demux type in a different column if needed
        if 'Sample_Type' in sample_row and sample_row['Sample_Type']:
            result['Demux_Sample_Type'] = sample_row['Sample_Type']
        
        # Process coverage data if matched
        if cov_match:
            pct_cov = cov_match.get('Exome_Percent_Coverage', None)
            
            # Debug: show what coverage value we got
            self.log(f"Processing coverage for {sample_id}: pct_cov = '{pct_cov}' (type: {type(pct_cov).__name__})")
            
            if pct_cov is not None and pct_cov != '':
                try:
                    pct_value = float(pct_cov)
                    result['Exome_Percent_Coverage'] = pct_value
                    
                    # Determine QC status
                    if pct_value >= self.threshold:
                        result['QC_Status'] = 'PASS'
                    else:
                        result['QC_Status'] = 'FAIL'
                    
                except (ValueError, TypeError) as e:
                    self.log(f"ERROR: Could not convert '{pct_cov}' to float: {e}")
                    result['Exome_Percent_Coverage'] = pct_cov
                    result['QC_Status'] = 'UNKNOWN'
            else:
                self.log(f"WARNING: Coverage value is None or empty for {sample_id}")
                result['Exome_Percent_Coverage'] = 'NA'
                result['QC_Status'] = 'NO_DATA'
        else:
            result['Exome_Percent_Coverage'] = 'NA'
            result['QC_Status'] = 'NOT_FOUND'
        
        return result
    
    def run(self):
        """Execute exome coverage QC check."""
        try:
            # Read demux file
            self.log(f"Reading demux file: {self.demux_file}")
            demux_samples = self.parse_demux()
            self.log(f"Loaded {len(demux_samples)} samples from demux")
            
            # Read coverage files (returns tumor_dict and normal_dict separately)
            self.log(f"Reading coverage files...")
            tumor_dict, normal_dict = self.read_coverage_files()
            
            if not tumor_dict and not normal_dict:
                print("ERROR: No coverage data found in tumor or normal files")
                return None
            
            # Process each sample in demux
            results = []
            controls_skipped = 0
            for sample_row in demux_samples:
                sample_id = sample_row.get('Sample_ID', '')
                
                # Skip control samples (NTC, NC, SC)
                if not self.is_sample(sample_id):
                    self.log(f"Skipping control sample: {sample_id}")
                    controls_skipped += 1
                    continue
                
                result = self.process_sample(sample_row, tumor_dict, normal_dict)
                results.append(result)
            
            if controls_skipped > 0:
                print(f"Skipped {controls_skipped} control samples (NTC, NC, SC)")
            
            # Save output using CSV writer
            self.save_results(results)
            
            print(f"Exome coverage QC saved: {self.output_file}")
            self.print_summary(results)
            
            return results
            
        except Exception as e:
            print(f"ERROR in ExomeCoverageQCChecker: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def save_results(self, results):
        """Save results to CSV file using csv.DictWriter."""
        if not results:
            print("WARNING: No results to save")
            return
        
        # Create output directory if needed
        self.output_file.parent.mkdir(parents=True, exist_ok=True)
        
        # Get all column names from first result
        fieldnames = list(results[0].keys())
        
        # Write CSV
        with open(self.output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
    
    def print_summary(self, results):
        """Print QC summary."""
        print("\n" + "="*70)
        print("EXOME COVERAGE QC SUMMARY")
        print(f"Coverage Threshold: >= {self.threshold}%")
        print("="*70)
        
        for row in results:
            sample = row['Sample_ID']
            specimen = row['Case_ID']
            exome_sample = row['Exome_Sample']
            coverage = row['Exome_Percent_Coverage']
            status = row['QC_Status']
            matched = row['Matched']
            
            if isinstance(coverage, (int, float)):
                print(f"{specimen:25s} | Coverage: {coverage:6.2f}% | {status:10s} | Matched: {matched}")
            else:
                print(f"{specimen:25s} | Coverage: {coverage:>10s} | {status:10s} | Matched: {matched}")
            
            if status == 'FAIL':
                print(f"    ⚠️  Coverage {coverage}% < {self.threshold}%")
            elif matched == 'No':
                print(f"    ℹ️  No coverage data found for this sample")
        
        print("="*70)
        
        # Overall stats
        pass_count = sum(1 for r in results if r['QC_Status'] == 'PASS')
        fail_count = sum(1 for r in results if r['QC_Status'] == 'FAIL')
        not_found_count = sum(1 for r in results if r['QC_Status'] == 'NOT_FOUND')
        no_data_count = sum(1 for r in results if r['QC_Status'] == 'NO_DATA')
        total = len(results)
        
        print(f"\nTotal: {total} | Pass: {pass_count} | Fail: {fail_count} | Not Found: {not_found_count} | No Data: {no_data_count}")


# Example usage / test code

class VCFOrganizer:
    """Organize VCF files into pass/pending folders based on QC results."""
    
    def __init__(self, vcf_dir, controls_qc_file, combined_qc_file, pass_dir=None, pending_dir=None):
        self.vcf_dir = Path(vcf_dir)
        self.controls_qc_file = Path(controls_qc_file)
        self.combined_qc_file = Path(combined_qc_file)
        
        # Set output directories
        # Default to vcf_dir parent if not specified
        default_output_dir = self.vcf_dir.parent
        self.pass_dir = Path(pass_dir) if pass_dir else default_output_dir / 'pass'
        self.pending_dir = Path(pending_dir) if pending_dir else default_output_dir / 'pending'

    
    def log(self, message):
        """Logging disabled."""
        pass
    
    def find_vcf_files(self, sample_id):
        """Find all VCF files where sample_id is contained in the filename."""
        vcf_files = []
        
        # Get all VCF files in directory
        all_vcfs = list(self.vcf_dir.glob('*.vcf*'))
        
        if not all_vcfs:
            self.log(f"No VCF files found in {self.vcf_dir}")
            return []
        
        # Simple substring matching
        for vcf_file in all_vcfs:
            if sample_id in vcf_file.name:
                vcf_files.append(vcf_file)
                self.log(f"  Match: {sample_id} in {vcf_file.name}")
        
        if not vcf_files:
            self.log(f"  No match: {sample_id} not found in any VCF filename")
        
        return vcf_files
    
    def copy_file(self, source_file, destination_dir):
        """Copy file to destination directory."""
        destination_dir.mkdir(parents=True, exist_ok=True)
        dest_file = destination_dir / source_file.name
        
        try:
            import shutil
            shutil.copy2(str(source_file), str(dest_file))
            self.log(f"Copied: {source_file.name} → {destination_dir.name}/")
            return True
        except Exception as e:
            print(f"ERROR copying {source_file.name}: {e}")
            return False
    
    def check_controls(self):
        """Check if all controls passed QC."""
        if not self.controls_qc_file.exists():
            print(f"WARNING: Controls QC file not found: {self.controls_qc_file}")
            return None
        
        controls_df = pd.read_csv(self.controls_qc_file)
        
        if 'Overall_QC' not in controls_df.columns:
            print("WARNING: Overall_QC column not found in controls_qc.csv")
            return None
        
        # Check if any control failed
        failed_controls = controls_df[controls_df['Overall_QC'] == 'FAIL']
        
        if len(failed_controls) > 0:
            print(f"\nControls QC: FAIL ({len(failed_controls)} controls failed)")
            for _, row in failed_controls.iterrows():
                print(f"  Failed control: {row['Sample_ID']}")
            return False
        else:
            print(f"\nControls QC: PASS (all {len(controls_df)} controls passed)")
            return True
    
    def organize_vcfs(self):
        """Organize VCF files based on QC results."""
        print("\n" + "="*70)
        print("VCF FILE ORGANIZATION")
        print("="*70)
        
        # Step 1: Check controls
        controls_pass = self.check_controls()
        
        if controls_pass is None:
            print("Cannot organize VCFs - controls QC file issue")
            return
        
        # Step 2: If any control failed, move ALL VCFs to pending
        if not controls_pass:
            print("\nCopying ALL VCF files to pending/ (controls failed)")
            all_vcfs = list(self.vcf_dir.glob('*.vcf*'))
            
            copied_count = 0
            for vcf_file in all_vcfs:
                if self.copy_file(vcf_file, self.pending_dir):
                    copied_count += 1
            
            print(f"\nMoved {copied_count} VCF files to pending/")
            print(f"Reason: One or more controls failed QC")
            return
        
        # Step 3: Controls passed - organize by Case_Level_QC
        if not self.combined_qc_file.exists():
            print(f"\nWARNING: Combined QC file not found: {self.combined_qc_file}")
            print("Cannot organize VCFs by case-level QC")
            return
        
        combined_df = pd.read_csv(self.combined_qc_file)
        
        if 'Case_Level_QC' not in combined_df.columns:
            print("WARNING: Case_Level_QC column not found in combined_qc.csv")
            return
        
        if 'Sample_ID' not in combined_df.columns:
            print("WARNING: Sample_ID column not found in combined_qc.csv")
            return
        
        print("\nOrganizing VCFs based on Case_Level_QC...")
        print("  Logic: If Sample_ID is in VCF filename, copy based on Case_Level_QC")
        print()
        
        # Track what we've copied
        copied_vcfs = set()  # Track which VCF files we've already copied
        pass_count = 0
        pending_count = 0
        
        # Process each sample in combined_qc.csv
        for _, row in combined_df.iterrows():
            sample_id = str(row['Sample_ID'])
            case_qc = row.get('Case_Level_QC', 'PENDING_REVIEW')
            
            self.log(f"Processing {sample_id}: Case_Level_QC = {case_qc}")
            
            # Determine destination folder
            if case_qc == 'PASS':
                dest_dir = self.pass_dir
            else:
                dest_dir = self.pending_dir
            
            # Find all VCF files where sample_id is contained in filename
            vcf_files = self.find_vcf_files(sample_id)
            
            # Copy each VCF file (skip if already copied)
            for vcf_file in vcf_files:
                if vcf_file not in copied_vcfs:
                    if self.copy_file(vcf_file, dest_dir):
                        copied_vcfs.add(vcf_file)
                        if case_qc == 'PASS':
                            pass_count += 1
                        else:
                            pending_count += 1
                else:
                    self.log(f"  Already copied: {vcf_file.name}")
        
        # Print summary
        print()
        print("="*70)
        print("VCF ORGANIZATION SUMMARY")
        print("="*70)
        print(f"VCF files copied to pass/:    {pass_count}")
        print(f"VCF files copied to pending/: {pending_count}")
        print(f"Total VCF files copied:       {pass_count + pending_count}")
        
        print(f"\nOutput directories:")
        print(f"  Pass:    {self.pass_dir}")
        print(f"  Pending: {self.pending_dir}")
        
        print("="*70)


def main():
    parser = argparse.ArgumentParser(
        description='Complete QC Pipeline - Validate controls, samples, HSMetrics, and Conpair'
    )
    
    # Required arguments
    parser.add_argument('-r', '--runqc', required=True, help='RunQC TSV file (controls)')
    parser.add_argument('-q', '--qc', required=True, help='QC TSV file (samples)')
    parser.add_argument('-d', '--demux', required=True, help='Demux CSV file')
    
    # Optional QC files
    parser.add_argument('-m', '--hsmetrics', help='HSMetrics summary CSV file (optional)')
    parser.add_argument('-c', '--conpair', help='Conpair summary CSV file (optional)')
    parser.add_argument('--exome-tumor', help='Exome tumor coverage percentages CSV file (optional)')
    parser.add_argument('--exome-normal', help='Exome normal coverage percentages CSV file (optional)')
    
    # Run ID
    parser.add_argument('--runid', help='Run ID to add to output files (optional)')
    
    # Thresholds
    parser.add_argument('-t', '--threshold', type=float, default=0.75, 
                       help='HSMetrics PCT_TARGET_BASES_250X threshold (default: 0.75)')
    parser.add_argument('--concordance-threshold', type=float, default=75.0,
                       help='Conpair concordance threshold percent (default: 75.0)')
    parser.add_argument('--exome-coverage-threshold', type=float, default=90.0,
                       help='Exome coverage percentage threshold (default: 90.0)')
    
    # Output files
    parser.add_argument('-oc', '--output-controls', default='controls_qc.csv', 
                       help='Controls output CSV (default: controls_qc.csv)')
    parser.add_argument('-os', '--output-samples', default='samples_qc.csv', 
                       help='Samples output CSV (default: samples_qc.csv)')
    parser.add_argument('-oh', '--output-hsmetrics', default='hsmetrics_qc.csv', 
                       help='HSMetrics output CSV (default: hsmetrics_qc.csv)')
    parser.add_argument('-op', '--output-conpair', default='conpair_qc.csv',
                       help='Conpair output CSV (default: conpair_qc.csv)')
    parser.add_argument('-oe', '--output-exome', default='exome_coverage_qc.csv',
                       help='Exome coverage output CSV (default: exome_coverage_qc.csv)')
    parser.add_argument('-om', '--output-merged', default='combined_qc.csv',
                       help='Merged output CSV (default: combined_qc.csv)')
    
    # Merge option
    parser.add_argument('--merge', action='store_true',
                       help='Merge samples, HSMetrics, and Conpair results into single file')
    
    # VCF organization
    parser.add_argument('--organize-vcfs', action='store_true',
                       help='Organize VCF files into pass/pending folders based on QC results')
    parser.add_argument('--vcf-dir', help='Directory containing VCF files to organize')
    parser.add_argument('--pass-dir', help='Directory for passing VCF files (default: <vcf-dir-parent>/pass)')
    parser.add_argument('--pending-dir', help='Directory for pending VCF files (default: <vcf-dir-parent>/pending)')
    
    # Verbose
    
    args = parser.parse_args()
    
    try:
        # Step 1: Check controls
        print("="*70)
        print("STEP 1: CONTROLS QC CHECK")
        print("="*70)
        
        controls_checker = ControlsQCChecker(
            runqc_file=args.runqc,
            demux_file=args.demux,
            output_file=args.output_controls,
            runid=args.runid,
            
        )
        
        controls_results = controls_checker.run()
        
        if controls_results is None:
            print("\nERROR: Controls QC check failed")
            sys.exit(1)
        
        all_passed = all(controls_results['Overall_QC'] == 'PASS')
        
        if not all_passed:
            print("\nControls QC FAILED - Stopping pipeline")
            print("One or more controls did not pass QC")
            sys.exit(1)
        
        print("\nAll controls PASSED - Proceeding to samples QC")
        
        # Step 2: Check samples
        print("\n" + "="*70)
        print("STEP 2: SAMPLES QC CHECK")
        print("="*70)
        
        samples_checker = SamplesQCChecker(
            qc_file=args.qc,
            demux_file=args.demux,
            output_file=args.output_samples,
            runid=args.runid,
            
        )
        
        samples_results = samples_checker.run()
        
        if samples_results is None:
            print("\nWARNING: No samples to check")
        
        # Step 3: Check HSMetrics (optional)
        if args.hsmetrics:
            print("\n" + "="*70)
            print("STEP 3: HSMETRICS QC CHECK")
            print("="*70)
            
            hsmetrics_checker = HSMetricsQCChecker(
                hsmetrics_file=args.hsmetrics,
                demux_file=args.demux,
                output_file=args.output_hsmetrics,
                threshold=args.threshold,
                runid=args.runid,
                
            )
            
            hsmetrics_results = hsmetrics_checker.run()
            
            if hsmetrics_results is None:
                print("\nWARNING: No samples to check in HSMetrics")
        
        # Step 4: Check Conpair (optional)
        if args.conpair:
            print("\n" + "="*70)
            print("STEP 4: CONPAIR QC CHECK")
            print("="*70)
            
            conpair_checker = ConpairQCChecker(
                conpair_file=args.conpair,
                demux_file=args.demux,
                output_file=args.output_conpair,
                concordance_threshold=args.concordance_threshold,
                runid=args.runid,
                
            )
            
            conpair_results = conpair_checker.run()
            
            if conpair_results is None:
                print("\nWARNING: No samples to check in Conpair")
        
        # Step 4.5: Check Exome Coverage (optional)
        if args.exome_tumor or args.exome_normal:
            print("\n" + "="*70)
            print("STEP 4.5: EXOME COVERAGE QC CHECK")
            print("="*70)
            
            # Both files should be provided together
            if not args.exome_tumor or not args.exome_normal:
                print("WARNING: Both --exome-tumor and --exome-normal should be provided")
                print("Skipping exome coverage QC")
            else:
                exome_checker = ExomeCoverageQCChecker(
                    tumor_file=args.exome_tumor,
                    normal_file=args.exome_normal,
                    demux_file=args.demux,
                    output_file=args.output_exome,
                    threshold=args.exome_coverage_threshold,
                    runid=args.runid,
                    
                )
                
                exome_results = exome_checker.run()
                
                if exome_results is None:
                    print("\nWARNING: No samples to check in Exome Coverage")
        
        # Step 5: Merge results (optional)
        if args.merge:
            print("\n" + "="*70)
            print("STEP 5: MERGING QC RESULTS")
            print("="*70)
            
            # Determine which files to merge
            merge_hsmetrics = args.hsmetrics is not None
            merge_conpair = args.conpair is not None
            
            # Read samples QC
            print(f"Reading {args.output_samples}...")
            merged_df = pd.read_csv(args.output_samples)
            print(f"  Base: {len(merged_df)} samples")
            
            # Debug: Show if critical columns exist
            if False:
                print(f"  Columns in samples_qc.csv: {list(merged_df.columns)}")
            if 'TM_Number' not in merged_df.columns:
                print(f"  WARNING: TM_Number column not found in {args.output_samples}")
                print(f"  Case_Level_QC will not be created")
            
            # Merge HSMetrics if it was run
            if merge_hsmetrics and Path(args.output_hsmetrics).exists():
                print(f"Merging {args.output_hsmetrics}...")
                hsmetrics_df = pd.read_csv(args.output_hsmetrics)
                
                # Select columns to merge
                hsmetrics_cols = ['Sample_ID', 'HSMetrics_Sample', 'PCT_TARGET_BASES_250X']
                if 'QC_Status' in hsmetrics_df.columns:
                    hsmetrics_df = hsmetrics_df.rename(columns={'QC_Status': 'HSMetrics_QC_Status'})
                    hsmetrics_cols.append('HSMetrics_QC_Status')
                if 'Matched' in hsmetrics_df.columns:
                    hsmetrics_cols.append('Matched')
                
                hsmetrics_merge = hsmetrics_df[[col for col in hsmetrics_cols if col in hsmetrics_df.columns]].copy()
                if 'Matched' in hsmetrics_merge.columns:
                    hsmetrics_merge = hsmetrics_merge.rename(columns={'Matched': 'HSMetrics_Matched'})
                
                merged_df = merged_df.merge(hsmetrics_merge, on='Sample_ID', how='left')
                print(f"  After HSMetrics: {len(merged_df)} rows, {len(merged_df.columns)} columns")
            
            # Merge Conpair if it was run
            if merge_conpair and Path(args.output_conpair).exists():
                print(f"Merging {args.output_conpair}...")
                conpair_df = pd.read_csv(args.output_conpair)
                
                # Select columns to merge
                conpair_cols = [
                    'Sample_ID',
                    'Conpair_Sample',
                    'Concordance_Percent',
                    'Concordance_Status'
                ]
                
                # Add optional columns with prefix
                if 'Matched' in conpair_df.columns:
                    conpair_cols.append('Matched')
                if 'Overall_QC' in conpair_df.columns:
                    conpair_cols.append('Overall_QC')
                if 'Failed_Metrics' in conpair_df.columns:
                    conpair_cols.append('Failed_Metrics')
                if 'Sample_Type' in conpair_df.columns:
                    conpair_cols.append('Sample_Type')
                
                conpair_merge = conpair_df[[col for col in conpair_cols if col in conpair_df.columns]].copy()
                
                # Rename to avoid conflicts
                rename_map = {
                    'Matched': 'Conpair_Matched',
                    'Overall_QC': 'Conpair_Overall_QC',
                    'Failed_Metrics': 'Conpair_Failed_Metrics',
                    'Sample_Type': 'Conpair_Sample_Type'
                }
                for old_name, new_name in rename_map.items():
                    if old_name in conpair_merge.columns:
                        conpair_merge = conpair_merge.rename(columns={old_name: new_name})
                
                merged_df = merged_df.merge(conpair_merge, on='Sample_ID', how='left')
                print(f"  After Conpair: {len(merged_df)} rows, {len(merged_df.columns)} columns")
            
            # Merge Exome Coverage if it was run
            merge_exome = (args.exome_tumor is not None and args.exome_normal is not None)
            if merge_exome and Path(args.output_exome).exists():
                print(f"Merging {args.output_exome}...")
                exome_df = pd.read_csv(args.output_exome)
                
                # Select columns to merge
                exome_cols = ['Sample_ID', 'Exome_Sample', 'Exome_Percent_Coverage']
                if 'Matched' in exome_df.columns:
                    exome_cols.append('Matched')
                if 'QC_Status' in exome_df.columns:
                    exome_cols.append('QC_Status')
                
                exome_merge = exome_df[[col for col in exome_cols if col in exome_df.columns]].copy()
                
                # Rename columns to avoid conflicts
                rename_map = {
                    'Matched': 'Exome_Matched',
                    'QC_Status': 'Exome_QC_Status'
                }
                for old_name, new_name in rename_map.items():
                    if old_name in exome_merge.columns:
                        exome_merge = exome_merge.rename(columns={old_name: new_name})
                
                merged_df = merged_df.merge(exome_merge, on='Sample_ID', how='left')
                print(f"  After Exome Coverage: {len(merged_df)} rows, {len(merged_df.columns)} columns")
            
            # Create Sample_Overall_QC and Case_Level_QC based on all QC statuses
            if 'TM_Number' in merged_df.columns and 'Sample_Level_QC' in merged_df.columns:
                print(f"\nCreating Sample_Overall_QC and Case_Level_QC...")
                
                # Determine which QC checks are available
                qc_checks = ['Sample_Level_QC']
                if 'HSMetrics_QC_Status' in merged_df.columns:
                    qc_checks.append('HSMetrics_QC_Status')
                if 'Conpair_Overall_QC' in merged_df.columns:
                    qc_checks.append('Conpair_Overall_QC')
                if 'Exome_QC_Status' in merged_df.columns:
                    qc_checks.append('Exome_QC_Status')
                
                print(f"  QC checks included: {', '.join(qc_checks)}")
                print(f"  All checks must PASS for Sample_Overall_QC = PASS")
                
                # Step 1: Create Sample_Overall_QC (per individual sample)
                def get_sample_overall_qc(row):
                    # Check if this sample passes all QC checks
                    sample_qc_pass = row.get('Sample_Level_QC') == 'PASS'
                    hsmetrics_qc_pass = row.get('HSMetrics_QC_Status', 'PASS') == 'PASS'
                    conpair_qc_pass = row.get('Conpair_Overall_QC', 'PASS') == 'PASS'
                    exome_qc_pass = row.get('Exome_QC_Status', 'PASS') == 'PASS'
                    
                    # All checks must pass for PASS, otherwise FAIL
                    if sample_qc_pass and hsmetrics_qc_pass and conpair_qc_pass and exome_qc_pass:
                        return 'PASS'
                    else:
                        return 'FAIL'
                
                # Apply to each row to get individual sample QC
                merged_df['Sample_Overall_QC'] = merged_df.apply(get_sample_overall_qc, axis=1)
                
                # Step 2: Create Case_Level_QC (per TM_Number)
                # Store TM_Number column separately to ensure it's preserved
                tm_number_col = merged_df['TM_Number'].copy()
                
                def update_case_qc(group):
                    # If ANY sample has Sample_Overall_QC = FAIL, case is PENDING_REVIEW
                    if (group['Sample_Overall_QC'] == 'FAIL').any():
                        group['Case_Level_QC'] = 'PENDING_REVIEW'
                    elif (group['Sample_Overall_QC'] == 'PASS').all():
                        group['Case_Level_QC'] = 'PASS'
                    else:
                        group['Case_Level_QC'] = 'NO_DATA'
                    return group
                
                # Apply to each TM_Number group - this preserves all columns including TM_Number
                merged_df = merged_df.groupby('TM_Number', as_index=False, group_keys=False).apply(update_case_qc)
                
                # Safety check: ensure TM_Number is still a column (not index)
                if 'TM_Number' not in merged_df.columns:
                    # If TM_Number became the index, reset it
                    if merged_df.index.name == 'TM_Number' or 'TM_Number' in str(merged_df.index.names):
                        merged_df = merged_df.reset_index()
                    else:
                        # Restore from our backup
                        merged_df['TM_Number'] = tm_number_col
                
                # Print summary
                sample_overall_summary = merged_df['Sample_Overall_QC'].value_counts()
                print(f"  Sample Overall QC summary:")
                for status, count in sample_overall_summary.items():
                    print(f"    {status}: {count} samples")
                
                # Only print case-level summary if TM_Number still exists
                if 'TM_Number' in merged_df.columns and 'Case_Level_QC' in merged_df.columns:
                    case_summary = merged_df.groupby('TM_Number')['Case_Level_QC'].first().value_counts()
                    print(f"  Case-level QC summary:")
                    for status, count in case_summary.items():
                        print(f"    {status}: {count} cases")
                else:
                    print(f"  Note: Case-level QC created but TM_Number column missing from output")
            else:
                # If no TM_Number column, still create Sample_Overall_QC without Case_Level_QC
                print(f"\nCreating Sample_Overall_QC (TM_Number column not found - skipping Case_Level_QC)...")
                
                # Determine which QC checks are available
                qc_checks = ['Sample_Level_QC']
                if 'HSMetrics_QC_Status' in merged_df.columns:
                    qc_checks.append('HSMetrics_QC_Status')
                if 'Conpair_Overall_QC' in merged_df.columns:
                    qc_checks.append('Conpair_Overall_QC')
                if 'Exome_QC_Status' in merged_df.columns:
                    qc_checks.append('Exome_QC_Status')
                
                if len(qc_checks) > 1:  # Only create if we have more than just Sample_Level_QC
                    print(f"  QC checks included: {', '.join(qc_checks)}")
                    print(f"  All checks must PASS for Sample_Overall_QC = PASS")
                    
                    # Create Sample_Overall_QC (per individual sample)
                    def get_sample_overall_qc(row):
                        # Check if this sample passes all QC checks
                        sample_qc_pass = row.get('Sample_Level_QC', 'PASS') == 'PASS'
                        hsmetrics_qc_pass = row.get('HSMetrics_QC_Status', 'PASS') == 'PASS'
                        conpair_qc_pass = row.get('Conpair_Overall_QC', 'PASS') == 'PASS'
                        exome_qc_pass = row.get('Exome_QC_Status', 'PASS') == 'PASS'
                        
                        # All checks must pass for PASS, otherwise FAIL
                        if sample_qc_pass and hsmetrics_qc_pass and conpair_qc_pass and exome_qc_pass:
                            return 'PASS'
                        else:
                            return 'FAIL'
                    
                    # Apply to each row to get individual sample QC
                    merged_df['Sample_Overall_QC'] = merged_df.apply(get_sample_overall_qc, axis=1)
                    
                    # Print summary
                    sample_overall_summary = merged_df['Sample_Overall_QC'].value_counts()
                    print(f"  Sample Overall QC summary:")
                    for status, count in sample_overall_summary.items():
                        print(f"    {status}: {count} samples")
            
            # Save merged file
            merged_output = Path(args.output_merged)
            merged_output.parent.mkdir(parents=True, exist_ok=True)
            merged_df.to_csv(merged_output, index=False)
            
            print(f"\nMerged QC report saved: {merged_output}")
            print(f"Total samples: {len(merged_df)}")
            print(f"Total columns: {len(merged_df.columns)}")
        
        # Step 6: Organize VCF files (if requested)
        if args.organize_vcfs:
            if not args.vcf_dir:
                print("\nERROR: --vcf-dir is required when using --organize-vcfs")
                sys.exit(1)
            
            if not args.merge:
                print("\nWARNING: --organize-vcfs requires --merge to be enabled")
                print("Enabling merge automatically...")
                # We need combined_qc.csv for organization
                if not Path(args.output_merged).exists():
                    print(f"ERROR: {args.output_merged} not found. Run with --merge first.")
                    sys.exit(1)
            
            print("\n" + "="*70)
            print("STEP 6: ORGANIZE VCF FILES")
            print("="*70)
            
            vcf_organizer = VCFOrganizer(
                vcf_dir=args.vcf_dir,
                controls_qc_file=args.output_controls,
                combined_qc_file=args.output_merged,
                pass_dir=args.pass_dir,
                pending_dir=args.pending_dir,
                
            )
            
            vcf_organizer.organize_vcfs()
        
        print("\n" + "="*70)
        print("QC PIPELINE COMPLETED")
        print("="*70)
        
    except Exception as e:
        print(f"ERROR: {str(e)}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()