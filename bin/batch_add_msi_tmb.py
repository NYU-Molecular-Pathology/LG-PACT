#!/usr/bin/env python3
"""
Process multiple VCF files to add both MSI and TMB values.

This script will:
1. Find all VCF files in a directory
2. Match them with MSI and TMB values from the TSV files
3. Add MSI and TMB headers to each VCF file
"""

import sys
import gzip
import csv
import argparse
from pathlib import Path


def read_annotation_values(annotation_file, annotation_type):
    """Read annotation values from TSV/CSV file into a dictionary."""
    annot_dict = {}
    
    # Determine delimiter
    with open(annotation_file, 'r') as f:
        first_line = f.readline().strip()
        delimiter = '\t' if '\t' in first_line else ','
    
    with open(annotation_file, 'r') as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        
        for row in reader:
            # Get sample ID from first column
            sample_id = None
            for key in row.keys():
                if 'sample' in key.lower() or key.strip() == list(row.keys())[0]:
                    sample_id = row[key].strip()
                    break
            
            if not sample_id:
                continue
            
            # Store all columns for this sample
            annot_dict[sample_id] = {}
            for key, value in row.items():
                if key and value:
                    column_name = key.strip()
                    # Special handling for MSI/TMB column -> rename to Score
                    if column_name == annotation_type:
                        column_name = 'Score'
                    annot_dict[sample_id][column_name] = value.strip()
    
    return annot_dict


def get_sample_name_from_vcf_path(vcf_path):
    """Extract sample name from VCF file path."""
    filename = Path(vcf_path).stem
    if filename.endswith('.vcf'):
        filename = filename[:-4]
    return filename


def is_gzipped(filename):
    """Check if file is gzipped."""
    return str(filename).endswith('.gz')


def open_file(filename, mode='r'):
    """Open file, handling gzipped files automatically."""
    if is_gzipped(filename):
        if 'b' not in mode:
            mode = mode.replace('r', 'rt').replace('w', 'wt')
        return gzip.open(filename, mode, compresslevel=6)
    else:
        return open(filename, mode)


def process_single_vcf(vcf_file, output_dir, msi_dict, tmb_dict, quiet=False):
    """Process a single VCF file."""
    try:
        vcf_filename = get_sample_name_from_vcf_path(vcf_file)
        
        result = {
            'file': vcf_file.name,
            'status': 'skipped',
            'msi_matched': False,
            'tmb_matched': False,
            'headers_added': 0,
            'message': ''
        }
        
        # Find MSI sample match
        msi_data = None
        if msi_dict:
            for sample_name in msi_dict.keys():
                if sample_name in vcf_filename:
                    msi_data = msi_dict[sample_name]
                    result['msi_matched'] = True
                    result['msi_sample'] = sample_name
                    break
        
        # Find TMB sample match
        tmb_data = None
        if tmb_dict:
            for sample_name in tmb_dict.keys():
                if sample_name in vcf_filename:
                    tmb_data = tmb_dict[sample_name]
                    result['tmb_matched'] = True
                    result['tmb_sample'] = sample_name
                    break
        
        # Check if we found at least one match
        if not msi_data and not tmb_data:
            result['message'] = 'No matching samples found'
            return result
        
        # Build headers
        all_headers = []
        
        if msi_data:
            for column, value in msi_data.items():
                if column.lower() not in ['sampleid', 'sample_id', 'sample', 'tumortype', 'tumor_type']:
                    all_headers.append(f"##MSI.{column}={value}\n")
        
        if tmb_data:
            for column, value in tmb_data.items():
                if column.lower() not in ['sampleid', 'sample_id', 'sample', 'tumortype', 'tumor_type', 'variantcaller', 'variant_caller']:
                    all_headers.append(f"##TMB.{column}={value}\n")
        
        result['headers_added'] = len(all_headers)
        
        # Write output file
        output_file = output_dir / vcf_file.name
        header_added = False
        
        with open_file(vcf_file, 'r') as infile, open_file(output_file, 'w') as outfile:
            for line in infile:
                if line.startswith('#CHROM') and not header_added:
                    for header in all_headers:
                        outfile.write(header)
                    header_added = True
                outfile.write(line)
        
        result['status'] = 'success'
        result['message'] = f"Added {len(all_headers)} header(s)"
        return result
        
    except Exception as e:
        return {
            'file': vcf_file.name,
            'status': 'error',
            'message': str(e),
            'msi_matched': False,
            'tmb_matched': False,
            'headers_added': 0
        }


def main():
    parser = argparse.ArgumentParser(description='Batch process VCF files to add MSI and/or TMB values')
    
    parser.add_argument('-d', '--vcf-dir', required=True, help='Directory containing VCF files')
    parser.add_argument('-m', '--msi', help='MSI TSV/CSV file with sample names and MSI values')
    parser.add_argument('-t', '--tmb', help='TMB TSV file with sample names and TMB values')
    parser.add_argument('-o', '--output-dir', required=True, help='Output directory for annotated VCF files')
    parser.add_argument('-q', '--quiet', action='store_true', help='Quiet mode - only show summary')
    
    args = parser.parse_args()
    
    # At least one of MSI or TMB must be provided
    if not args.msi and not args.tmb:
        print("Error: At least one of --msi or --tmb must be provided")
        sys.exit(1)
    
    vcf_dir = Path(args.vcf_dir)
    output_dir = Path(args.output_dir)
    
    # Validate input directory exists
    if not vcf_dir.exists():
        print(f"Error: VCF directory not found: {vcf_dir}")
        sys.exit(1)
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Read MSI values if provided
    msi_dict = None
    if args.msi:
        if not Path(args.msi).exists():
            print(f"Error: MSI file not found: {args.msi}")
            sys.exit(1)
        print(f"Reading MSI values from: {args.msi}")
        msi_dict = read_annotation_values(args.msi, 'MSI')
        print(f"Found {len(msi_dict)} samples with MSI values")
    
    # Read TMB values if provided
    tmb_dict = None
    if args.tmb:
        if not Path(args.tmb).exists():
            print(f"Error: TMB file not found: {args.tmb}")
            sys.exit(1)
        print(f"Reading TMB values from: {args.tmb}")
        tmb_dict = read_annotation_values(args.tmb, 'TMB')
        print(f"Found {len(tmb_dict)} samples with TMB values")
    
    print()
    
    # Find all VCF files
    vcf_files = list(vcf_dir.glob("*.vcf")) + list(vcf_dir.glob("*.vcf.gz"))
    
    if not vcf_files:
        print(f"No VCF files found in {vcf_dir}")
        sys.exit(1)
    
    print(f"Found {len(vcf_files)} VCF file(s)")
    print(f"{'='*60}\n")
    
    # Process VCF files
    results = []
    for vcf_file in vcf_files:
        result = process_single_vcf(vcf_file, output_dir, msi_dict, tmb_dict, args.quiet)
        results.append(result)
    
    # Analyze results
    success = sum(1 for r in results if r['status'] == 'success')
    skipped = sum(1 for r in results if r['status'] == 'skipped')
    errors = sum(1 for r in results if r['status'] == 'error')
    partial_matches = sum(1 for r in results if r['status'] == 'success' and (
        (msi_dict and not r['msi_matched']) or (tmb_dict and not r['tmb_matched'])
    ))
    
    # Print individual results if not in quiet mode
    if not args.quiet:
        for result in results:
            status_symbol = '✓' if result['status'] == 'success' else '✗'
            print(f"{status_symbol} {result['file']}: {result['message']}")
            if result['status'] == 'success':
                if result.get('msi_matched'):
                    print(f"  MSI: {result.get('msi_sample')}")
                if result.get('tmb_matched'):
                    print(f"  TMB: {result.get('tmb_sample')}")
        print()
    
    # Print summary
    print(f"{'='*60}")
    print("Summary:")
    print(f"  Total:      {len(vcf_files)}")
    print(f"  Success:    {success}")
    if partial_matches > 0:
        print(f"  Partial:    {partial_matches}")
    if skipped > 0:
        print(f"  Skipped:    {skipped}")
    if errors > 0:
        print(f"  Errors:     {errors}")
    print(f"  Output:     {output_dir}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()