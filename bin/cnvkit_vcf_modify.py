#!/usr/bin/env python3
"""
Modify SVTYPE in CNVkit VCF file based on CN and FOLD_CHANGE_LOG values:
- If CN present and CN > 5: SVTYPE = Amplification
- If CN present and 3 <= CN <= 5: SVTYPE = Gain
- If CN present and CN = 2: SVTYPE = Neutral
- If FOLD_CHANGE_LOG <= -1.1: SVTYPE = Loss and add CN=0
- If -1.1 < FOLD_CHANGE_LOG < 0: SVTYPE = LOH and add CN=1

Optionally change sample header name.
"""

import sys
import re
import argparse

def extract_fold_change_log(info_field):
    """Extract FOLD_CHANGE_LOG value from INFO field"""
    match = re.search(r'FOLD_CHANGE_LOG=([-\d.]+)', info_field)
    if match:
        return float(match.group(1))
    return None

def modify_vcf(input_file, output_file, sample_name=None):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Handle header lines
            if line.startswith('#'):
                # Change sample name if requested
                if sample_name and line.startswith('#CHROM'):
                    line = line.rstrip()
                    # Replace the last column header (sample name)
                    parts = line.split('\t')
                    if len(parts) >= 10:
                        parts[-1] = sample_name
                        line = '\t'.join(parts) + '\n'
                    else:
                        line = line + '\n'
                outfile.write(line)
                continue
            
            # Process variant lines
            fields = line.strip().split('\t')
            if len(fields) < 10:
                outfile.write(line)
                continue
            
            # Extract FORMAT and sample fields
            format_field = fields[8]
            sample_field = fields[9]
            info_field = fields[7]
            
            format_tags = format_field.split(':')
            sample_values = sample_field.split(':')
            
            # Check if CN is in FORMAT
            cn_value = None
            has_cn = False
            if 'CN' in format_tags:
                has_cn = True
                cn_index = format_tags.index('CN')
                if cn_index < len(sample_values):
                    try:
                        cn_value = int(sample_values[cn_index])
                    except (ValueError, IndexError):
                        cn_value = None
            
            # Extract FOLD_CHANGE_LOG
            fold_change_log = extract_fold_change_log(info_field)
            
            # Determine new SVTYPE and whether to add CN
            new_svtype = None
            add_cn = False
            cn_to_add = None
            
            if has_cn and cn_value is not None:
                # Rules for when CN is present
                if cn_value > 5:
                    new_svtype = 'Amplification'
                elif 3 <= cn_value <= 5:
                    new_svtype = 'Gain'
                elif cn_value == 2:
                    new_svtype = 'Neutral'
                # If CN < 2, keep original SVTYPE
            else:
                # Rules for when CN is missing - based on FOLD_CHANGE_LOG
                if fold_change_log is not None:
                    if fold_change_log <= -1.1:
                        new_svtype = 'Loss'
                        add_cn = True
                        cn_to_add = 0
                    elif -1.1 < fold_change_log < 0:
                        new_svtype = 'LOH'
                        add_cn = True
                        cn_to_add = 1
            
            # If no changes needed, write original line
            if new_svtype is None and not add_cn:
                outfile.write(line)
                continue
            
            # Modify INFO field SVTYPE if needed
            if new_svtype:
                info_field = re.sub(r'SVTYPE=[^;]+', f'SVTYPE={new_svtype}', info_field)
                fields[7] = info_field
            
            # Add CN to FORMAT and sample if needed
            if add_cn and cn_to_add is not None:
                # Add CN to FORMAT
                if 'CN' not in format_tags:
                    format_tags.append('CN')
                    sample_values.append(str(cn_to_add))
                    fields[8] = ':'.join(format_tags)
                    fields[9] = ':'.join(sample_values)
            
            # Write modified line
            outfile.write('\t'.join(fields) + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Modify SVTYPE in CNVkit VCF based on CN and FOLD_CHANGE_LOG values'
    )
    parser.add_argument('--input', '-i', help='Input VCF file')
    parser.add_argument('--output', '-o', help='Output VCF file')
    parser.add_argument('--sample-header', '-s', help='New sample name for the header (optional)')
    
    args = parser.parse_args()
    
    modify_vcf(args.input, args.output, args.sample_header)
    
    print(f"Modified VCF written to: {args.output}")
    if args.sample_header:
        print(f"Sample header changed to: {args.sample_header}")