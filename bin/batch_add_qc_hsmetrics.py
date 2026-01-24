#!/usr/bin/env python3
"""
Batch script to add QC and HSMetrics headers to VCF files.

Parses QC TSV file, HSMetrics HTML file, and tumor-normal pairing CSV,
then adds comprehensive headers to matching VCF files.
Supports both .vcf and .vcf.gz files, outputs .vcf.gz files.
"""

import argparse
import csv
import glob
import gzip
import os
import re
import sys
from html.parser import HTMLParser


# Constants
CONTROL_SAMPLES = {'NC', 'NTC', 'SC'}
SAMPLE_ID_PATTERN = r'([A-Za-z0-9]+-\d+(?:-[A-Za-z0-9]+)?_B\d+-\d+)'
BASE_ID_PATTERN = r'([A-Za-z0-9]+-\d+)'

HSMETRICS_FIELDS = [
    'ON_BAIT_BASES',
    'MEAN_BAIT_COVERAGE',
    'PCT_TARGET_BASES_50X',
    'PCT_TARGET_BASES_250X',
    'PCT_USABLE_BASES_ON_BAIT',
    'PCT_USABLE_BASES_ON_TARGET',
    'TOTAL_READS',
    'PF_READS',
    'ON_TARGET_BASES',
    'MEAN_TARGET_COVERAGE',
    'FOLD_80_BASE_PENALTY',
    'AT_DROPOUT',
    'GC_DROPOUT'
]


class HSMetricsParser(HTMLParser):
    """Parser for HSMetrics HTML report files."""
    
    def __init__(self):
        super().__init__()
        self.in_table = False
        self.in_row = False
        self.in_header = False
        self.current_row = []
        self.headers = []
        self.data = []
        self.current_cell = []
        
    def handle_starttag(self, tag, attrs):
        if tag == 'table':
            self.in_table = True
        elif tag == 'tr' and self.in_table:
            self.in_row = True
            self.current_row = []
        elif tag == 'th' and self.in_row:
            self.in_header = True
            self.current_cell = []
        elif tag == 'td' and self.in_row:
            self.current_cell = []
            
    def handle_endtag(self, tag):
        if tag == 'table':
            self.in_table = False
        elif tag == 'tr':
            if self.in_row and self.current_row:
                if self.headers:
                    self.data.append(self.current_row)
                else:
                    self.headers = self.current_row
            self.in_row = False
        elif tag in ('th', 'td'):
            if self.current_cell:
                self.current_row.append(''.join(self.current_cell).strip())
            self.in_header = tag == 'th' and False
                
    def handle_data(self, data):
        if self.in_row:
            self.current_cell.append(data)


def extract_sample_id(full_name):
    """Extract sample ID from the full sample name."""
    match = re.search(SAMPLE_ID_PATTERN, full_name)
    if match:
        return match.group(1)
    match = re.search(r'([A-Za-z0-9]+-\d+(?:-[A-Za-z0-9]+)?)', full_name)
    return match.group(1) if match else full_name


def extract_base_id(sample_id):
    """Extract base ID from sample ID (without block number)."""
    match = re.match(BASE_ID_PATTERN, sample_id)
    return match.group(1) if match else sample_id


def parse_qc_value(value_str):
    """Extract numeric value from QC field."""
    return value_str.strip().split()[0] if value_str.strip() else ""


def parse_percentage(value_str):
    """Extract percentage value from string like '0.9% Pass'."""
    match = re.match(r'([\d.]+)%', value_str.strip())
    return str(float(match.group(1))) if match else ""


def parse_coverage(coverage_str):
    """Parse mean and median coverage from 'mean,median Pass' format."""
    match = re.match(r'([\d.]+),([\d.]+)', coverage_str)
    return (match.group(1), match.group(2)) if match else ("", "")


def read_file_safely(filepath, file_type):
    """Check if file exists and raise informative error if not."""
    if not os.path.exists(filepath):
        print(f"Error: {file_type} file not found: {filepath}", file=sys.stderr)
        sys.exit(1)


def read_tumor_normal_pairs(filepath):
    """Read tumor-normal pairing from CSV file."""
    read_file_safely(filepath, "Pairing CSV")
    pairs = {}
    
    try:
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                tumor_id = extract_sample_id(row['Tumor'])
                normal_id = extract_sample_id(row['Normal'])
                base_id = extract_base_id(tumor_id)
                pairs[base_id] = {'tumor': tumor_id, 'normal': normal_id}
    except KeyError as e:
        print(f"Error: Missing required column in {filepath}: {e}", file=sys.stderr)
        sys.exit(1)
    
    return pairs


def read_qc_data(filepath):
    """Read QC metrics from TSV file."""
    read_file_safely(filepath, "QC TSV")
    qc_data = {}
    
    try:
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                qc_data[row['Sample']] = row
    except KeyError as e:
        print(f"Error: Missing required column in {filepath}: {e}", file=sys.stderr)
        sys.exit(1)
    
    return qc_data


def read_hsmetrics_data(filepath):
    """Read HSMetrics from HTML file."""
    read_file_safely(filepath, "HSMetrics HTML")
    hsmetrics_data = {}
    
    with open(filepath, 'r') as f:
        html_content = f.read()
    
    parser = HSMetricsParser()
    parser.feed(html_content)
    
    for row in parser.data:
        if row and row[0] not in CONTROL_SAMPLES:
            sample_id = row[0]
            metrics = {
                header: row[i] 
                for i, header in enumerate(parser.headers[1:], 1) 
                if i < len(row)
            }
            hsmetrics_data[sample_id] = metrics
    
    return hsmetrics_data


def generate_metric_headers(sample_type, sample_id, qc_data, prefix="FLAGSTAT"):
    """Generate FLAGSTAT headers for a single sample (tumor or normal)."""
    if sample_id not in qc_data:
        return []
    
    qc = qc_data[sample_id]
    mean_cov, median_cov = parse_coverage(qc['Average and Median Coverage'])
    
    return [
        f"##{prefix}.{sample_type}.MAPPEDREADS={parse_qc_value(qc['Mapped Reads QC'])}",
        f"##{prefix}.{sample_type}.DEDUPLICATEDREADS={parse_qc_value(qc['Deduplicated Reads QC'])}",
        f"##{prefix}.{sample_type}.MEANCOVERAGE={mean_cov}",
        f"##{prefix}.{sample_type}.MEDIANCOVERAGE={median_cov}",
        f"##{prefix}.{sample_type}.PCT_TARGET_50X_COVERAGE={parse_percentage(qc['Targets with <50 Coverage'])}"
    ]


def generate_flagstat_headers(tumor_id, normal_id, qc_data):
    """Generate FLAGSTAT VCF headers for a tumor-normal pair."""
    return (generate_metric_headers('TUMOR', tumor_id, qc_data) +
            generate_metric_headers('NORMAL', normal_id, qc_data))


def find_hsmetrics_sample(base_id, hsmetrics_data, has_block):
    """Find HSMetrics sample matching criteria."""
    for sample_id, metrics in hsmetrics_data.items():
        if extract_base_id(sample_id) == base_id:
            if has_block and len(sample_id) > len(base_id):
                return metrics
            elif not has_block and sample_id == base_id:
                return metrics
    return None


def generate_hsmetrics_sample_headers(sample_type, metrics):
    """Generate HSMetrics headers for a single sample."""
    if not metrics:
        return []
    
    headers = []
    for field in HSMETRICS_FIELDS:
        value = metrics.get(field, '')
        if field == 'PCT_TARGET_BASES_50X' and not value:
            value = metrics.get('PCT_TARGET_BASES_250X', '')
        headers.append(f"##HSMETRICS.{sample_type}.{field}={value}")
    
    return headers


def generate_hsmetrics_headers(tumor_id, normal_id, hsmetrics_data):
    """Generate HSMETRICS VCF headers for a tumor-normal pair."""
    tumor_base = extract_base_id(tumor_id)
    normal_base = extract_base_id(normal_id)
    
    tumor_metrics = find_hsmetrics_sample(tumor_base, hsmetrics_data, has_block=True)
    normal_metrics = find_hsmetrics_sample(normal_base, hsmetrics_data, has_block=False)
    
    return (generate_hsmetrics_sample_headers('TUMOR', tumor_metrics) +
            generate_hsmetrics_sample_headers('NORMAL', normal_metrics))


def find_insertion_index(lines, insert_after):
    """Find the index where headers should be inserted."""
    if insert_after == 'msi_tmb':
        for i in range(len(lines) - 1, -1, -1):
            if lines[i].startswith('##MSI') or lines[i].startswith('##TMB'):
                return i + 1
        insert_after = 'fileformat'
    
    if insert_after == 'contig':
        for i in range(len(lines) - 1, -1, -1):
            if lines[i].startswith('##contig'):
                return i + 1
        insert_after = 'fileformat'
    
    for i, line in enumerate(lines):
        if line.startswith('##fileformat'):
            return i + 1
    
    return 1


def read_vcf_file(vcf_path):
    """Read VCF file (handles both .vcf and .vcf.gz)."""
    if vcf_path.endswith('.gz'):
        with gzip.open(vcf_path, 'rt') as f:
            return f.readlines()
    else:
        with open(vcf_path, 'r') as f:
            return f.readlines()


def write_vcf_file(output_path, lines):
    """Write VCF file as gzipped."""
    if not output_path.endswith('.gz'):
        output_path = output_path + '.gz'
    
    with gzip.open(output_path, 'wt') as f:
        f.writelines(lines)
    
    return output_path


def add_headers_to_vcf(vcf_path, headers, output_path, insert_after='msi_tmb'):
    """Add headers to a VCF file."""
    lines = read_vcf_file(vcf_path)
    
    insert_index = find_insertion_index(lines, insert_after)
    header_lines = [f"{h}\n" for h in headers]
    new_lines = lines[:insert_index] + header_lines + lines[insert_index:]
    
    final_output_path = write_vcf_file(output_path, new_lines)
    return final_output_path


def get_base_name(filepath):
    """Get base name without .vcf or .vcf.gz extension."""
    basename = os.path.basename(filepath)
    if basename.endswith('.vcf.gz'):
        return basename[:-7]
    elif basename.endswith('.vcf'):
        return basename[:-4]
    return basename


def find_matching_vcf_files(vcf_dir, base_id, pattern):
    """Find VCF files matching the base ID (both .vcf and .vcf.gz)."""
    all_vcfs = []
    
    # Search for both .vcf and .vcf.gz files
    for ext in ['*.vcf', '*.vcf.gz']:
        all_vcfs.extend(glob.glob(os.path.join(vcf_dir, ext)))
    
    # Filter by base_id in filename
    return [vcf for vcf in all_vcfs if base_id in os.path.basename(vcf)]


def process_vcf_files(vcf_files, headers, args):
    """Process a list of VCF files with given headers."""
    processed = 0
    errors = 0
    
    for vcf_path in vcf_files:
        vcf_basename = os.path.basename(vcf_path)
        base_name = get_base_name(vcf_path)
        
        if args.output_dir:
            output_path = os.path.join(args.output_dir, f"{base_name}{args.suffix}.vcf.gz")
        else:
            output_path = vcf_path.replace('.vcf.gz', f'{args.suffix}.vcf.gz').replace('.vcf', f'{args.suffix}.vcf.gz')
        
        if args.dry_run:
            print(f"[DRY RUN] {vcf_basename} -> {os.path.basename(output_path)} ({len(headers)} headers)")
        else:
            try:
                final_path = add_headers_to_vcf(vcf_path, headers, output_path, args.insert_after)
                print(f"✓ {vcf_basename} -> {os.path.basename(final_path)}")
                processed += 1
            except Exception as e:
                print(f"✗ {vcf_basename}: {e}", file=sys.stderr)
                errors += 1
    
    return processed, errors


def main():
    """Main function to batch process VCF files."""
    parser = argparse.ArgumentParser(
        description='Batch add FLAGSTAT and HSMETRICS headers to VCF files.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -q qc.tsv -m hsmetrics.html -p samples.csv -d vcf_files/
  %(prog)s -q qc.tsv -m hsmetrics.html -p samples.csv -d vcfs/ -o output/

Input:
  Accepts both .vcf and .vcf.gz files
  
Output:
  Adds 36 headers per VCF file (10 FLAGSTAT + 26 HSMETRICS)
  All output files are gzip compressed (.vcf.gz)
        """
    )
    
    parser.add_argument('-q', '--qc-file', required=True, help='Path to QC TSV file')
    parser.add_argument('-m', '--hsmetrics-file', required=True, help='Path to HSMetrics HTML file')
    parser.add_argument('-p', '--pairing-file', required=True, help='Path to tumor-normal pairing CSV')
    parser.add_argument('-d', '--vcf-dir', required=True, help='Directory containing VCF files (.vcf or .vcf.gz)')
    parser.add_argument('-o', '--output-dir', help='Output directory for processed VCFs')
    parser.add_argument('--pattern', default='*.vcf*', help='Glob pattern for VCF files (default: *.vcf*)')
    parser.add_argument('--insert-after', choices=['fileformat', 'msi_tmb', 'contig'], 
                       default='msi_tmb', help='Where to insert headers (default: msi_tmb)')
    parser.add_argument('--suffix', default='_qc', help='Suffix for output files (default: _qc)')
    parser.add_argument('--dry-run', action='store_true', help='Preview without modifying files')
    
    args = parser.parse_args()
    
    if args.output_dir:
        os.makedirs(args.output_dir, exist_ok=True)
    
    # Read all input files
    qc_data = read_qc_data(args.qc_file)
    hsmetrics_data = read_hsmetrics_data(args.hsmetrics_file)
    pairs = read_tumor_normal_pairs(args.pairing_file)
    
    # Process each pair
    total_processed = 0
    total_skipped = 0
    total_errors = 0
    
    for base_id, pair_info in pairs.items():
        vcf_files = find_matching_vcf_files(args.vcf_dir, base_id, args.pattern)
        
        if not vcf_files:
            print(f"⚠  No VCF files found for {base_id}", file=sys.stderr)
            total_skipped += 1
            continue
        
        # Generate all headers
        flagstat_headers = generate_flagstat_headers(
            pair_info['tumor'], pair_info['normal'], qc_data
        )
        hsmetrics_headers = generate_hsmetrics_headers(
            pair_info['tumor'], pair_info['normal'], hsmetrics_data
        )
        all_headers = flagstat_headers + hsmetrics_headers
        
        if not all_headers:
            print(f"⚠  No headers generated for {base_id}", file=sys.stderr)
            total_errors += 1
            continue
        
        # Process VCF files
        processed, errors = process_vcf_files(vcf_files, all_headers, args)
        total_processed += processed
        total_errors += errors
    
    # Print summary
    print(f"\n{'='*50}")
    print(f"{'DRY RUN ' if args.dry_run else ''}SUMMARY")
    print(f"{'='*50}")
    print(f"Pairs:      {len(pairs)}")
    print(f"Processed:  {total_processed}")
    print(f"Skipped:    {total_skipped}")
    print(f"Errors:     {total_errors}")
    print(f"{'='*50}")


if __name__ == '__main__':
    main()