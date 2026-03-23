#!/usr/bin/env python3
"""
Combined QC Parser and MySQL Upload (Direct Connection - No SSH)
Parses paired tumor/normal QC CSV and uploads to MySQL database via direct connection
"""
import csv
import pymysql
import pymysql.cursors
import sys
import argparse
import os

# =====================================================
# CONFIGURATION
# =====================================================
try:
    from db_config import DB_CONFIG
    try:
        from db_config import COMBINED_TABLE_NAME as TABLE_NAME
    except ImportError:
        from db_config import TABLE_NAME
    print("✓ Using configuration from config.py")
except ImportError:
    print("⚠ config.py not found, using default configuration")
    DB_CONFIG = {
        'host': 'localhost',
        'user': 'your_username',
        'password': 'your_password',
        'database': 'qc_database',
        'port': 3306
    }
    TABLE_NAME = 'sampleqc_metrics'


# =====================================================
# PARSING FUNCTIONS
# =====================================================

def parse_combined_qc_csv(input_file):
    """
    Parse paired tumor/normal QC CSV file
    Returns list of combined records (one per case)
    """
    parsed_data = []
    
    try:
        with open(input_file, 'r') as infile:
            reader = csv.DictReader(infile)
            rows = list(reader)
            
            # Group rows by TM_number (case_id)
            cases = {}
            for row in rows:
                tm_number = row.get('TM_Number', '').strip()
                sample_type = row.get('Sample_Type', '').strip()
                
                if tm_number not in cases:
                    cases[tm_number] = {}
                
                cases[tm_number][sample_type] = row
            
            # Process each case
            for tm_number, samples in cases.items():
                tumor_row  = samples.get('Tumor', {})
                normal_row = samples.get('Normal', {})
                
                if not tumor_row and not normal_row:
                    continue
                
                def get_value(row, key):
                    value = row.get(key, '').strip()
                    return 'NA' if value == '' else value
                
                def clean_percentage(value):
                    if value == 'NA' or value == '':
                        return 'NA'
                    return value.replace('%', '').strip()
                
                def clean_snp_overlap(value):
                    if value == 'NA' or value == '':
                        return 'NA'
                    if '%' in value:
                        return value.split('%')[0].strip()
                    return value.strip()
                
                combined_record = {
                    'run_id':                      get_value(tumor_row or normal_row, 'Run_ID'),
                    'case_id':                     tm_number if tm_number else 'NA',
                    'tumor_id':                    get_value(tumor_row, 'Case_ID'),
                    'normal_id':                   get_value(normal_row, 'Case_ID'),
                    'status':                      get_value(tumor_row or normal_row, 'Case_Level_QC'),
                    'tumor_status':                get_value(tumor_row, 'Sample_Overall_QC'),
                    'normal_status':               get_value(normal_row, 'Sample_Overall_QC'),
                    'tumor_mapped_reads':          get_value(tumor_row, 'Mapped Reads QC_Value'),
                    'tumor_dedup_reads':           get_value(tumor_row, 'Deduplicated Reads QC_Value'),
                    'tumor_pct_target_ge_50x':     clean_percentage(get_value(tumor_row, 'Targets with <50 Coverage_Value')),
                    'tumor_mean_coverage':         get_value(tumor_row, 'Average_Coverage'),
                    'tumor_median_coverage':       get_value(tumor_row, 'Median_Coverage'),
                    'tumor_snp_overlap':           clean_snp_overlap(get_value(tumor_row, 'Overlapped HOMO SNPs_Value')),
                    'tumor_concordance':           get_value(tumor_row, 'Concordance_Percent'),
                    'normal_mapped_reads':         get_value(normal_row, 'Mapped Reads QC_Value'),
                    'normal_dedup_reads':          get_value(normal_row, 'Deduplicated Reads QC_Value'),
                    'normal_pct_target_ge_50x':    clean_percentage(get_value(normal_row, 'Targets with <50 Coverage_Value')),
                    'normal_mean_coverage':        get_value(normal_row, 'Average_Coverage'),
                    'normal_median_coverage':      get_value(normal_row, 'Median_Coverage'),
                    'normal_snp_overlap':          clean_snp_overlap(get_value(normal_row, 'Overlapped HOMO SNPs_Value')),
                    'normal_concordance':          get_value(tumor_row, 'Concordance_Percent'),
                    'tumor_pct_target_bases_250x': clean_percentage(get_value(tumor_row, 'PCT_TARGET_BASES_250X')),
                    'tumor_exome_coverage':        get_value(tumor_row, 'Exome_Percent_Coverage'),
                    'normal_pct_target_bases_250x': clean_percentage(get_value(normal_row, 'PCT_TARGET_BASES_250X')),
                    'normal_exome_coverage':       get_value(normal_row, 'Exome_Percent_Coverage'),
                }
                
                parsed_data.append(combined_record)
        
        print(f"✓ Parsed {len(parsed_data)} case(s) from {input_file}")
        return parsed_data
        
    except FileNotFoundError:
        print(f"✗ Error: File '{input_file}' not found", file=sys.stderr)
        sys.exit(1)
    except KeyError as e:
        print(f"✗ Error: Required column not found: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"✗ Error parsing file: {e}", file=sys.stderr)
        sys.exit(1)


def save_parsed_csv(parsed_data, input_file):
    """Save parsed data to CSV file for reference"""
    try:
        base_name = os.path.basename(input_file)
        name_without_ext = os.path.splitext(base_name)[0]
        output_file = f'{name_without_ext}_parsed.csv'
        
        if parsed_data:
            fieldnames = [
                'run_id', 'case_id', 'tumor_id', 'normal_id', 'status',
                'tumor_status', 'normal_status',
                'tumor_mapped_reads', 'tumor_dedup_reads', 'tumor_pct_target_ge_50x',
                'tumor_mean_coverage', 'tumor_median_coverage', 'tumor_snp_overlap',
                'tumor_concordance','normal_mapped_reads', 'normal_dedup_reads', 'normal_pct_target_ge_50x',
                'normal_mean_coverage', 'normal_median_coverage', 'normal_snp_overlap',
                'normal_concordance','tumor_pct_target_bases_250x', 'tumor_exome_coverage',
                'normal_pct_target_bases_250x', 'normal_exome_coverage'
            ]
            
            with open(output_file, 'w', newline='') as outfile:
                writer = csv.DictWriter(outfile, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(parsed_data)
            
            print(f"✓ Saved parsed CSV to: {output_file}")
            return output_file
        
    except Exception as e:
        print(f"⚠ Warning: Could not save parsed CSV: {e}")
        return None


# =====================================================
# DATABASE FUNCTIONS
# =====================================================

def create_connection(config):
    """Create a database connection using PyMySQL"""
    try:
        connection = pymysql.connect(
            host=config['host'],
            user=config['user'],
            password=config['password'],
            database=config['database'],
            port=config.get('port', 3306),
            cursorclass=pymysql.cursors.DictCursor
        )
        print(f"✓ Connected to MySQL Server")
        print(f"✓ Using database: {config['database']}")
        return connection

    except pymysql.MySQLError as e:
        print(f"✗ Error connecting to MySQL: {e}", file=sys.stderr)
        sys.exit(1)


def verify_table_exists(connection, table_name):
    """Verify table exists - exit if it doesn't"""
    try:
        with connection.cursor() as cursor:
            cursor.execute(f"SHOW TABLES LIKE '{table_name}'")
            result = cursor.fetchone()
        
        if not result:
            print(f"✗ Error: Table '{table_name}' does not exist", file=sys.stderr)
            print("This script only works with existing tables.")
            sys.exit(1)
        
        print(f"✓ Table '{table_name}' exists")

    except pymysql.MySQLError as e:
        print(f"✗ Error checking table: {e}", file=sys.stderr)
        sys.exit(1)


def check_for_duplicates(connection, table_name, parsed_data):
    """Check if data already exists in table"""
    try:
        case_ids = set(row['case_id'] for row in parsed_data if row['case_id'] != 'NA')
        duplicates = []

        with connection.cursor() as cursor:
            for case_id in case_ids:
                cursor.execute(
                    f"SELECT COUNT(*) as count FROM {table_name} WHERE case_id = %s",
                    (case_id,)
                )
                count = cursor.fetchone()['count']
                if count > 0:
                    duplicates.append((case_id, count))

        return duplicates

    except Exception as e:
        print(f"⚠ Could not check for duplicates: {e}")
        return []


def delete_existing_cases(connection, table_name, case_ids):
    """Delete existing data for specific case_ids"""
    try:
        placeholders = ','.join(['%s'] * len(case_ids))
        query = f"DELETE FROM {table_name} WHERE case_id IN ({placeholders})"

        with connection.cursor() as cursor:
            cursor.execute(query, list(case_ids))
            deleted = cursor.rowcount

        connection.commit()
        print(f"✓ Deleted {deleted} existing records for {len(case_ids)} case_id(s)")

    except pymysql.MySQLError as e:
        connection.rollback()
        print(f"✗ Error deleting data: {e}", file=sys.stderr)
        sys.exit(1)


def insert_parsed_data(connection, parsed_data, table_name):
    """Insert parsed data into the database"""
    insert_query = f"""
    INSERT INTO {table_name} 
    (run_id, case_id, tumor_id, normal_id, status,
     tumor_status, normal_status,
     tumor_mapped_reads, tumor_dedup_reads, tumor_pct_target_ge_50x,
     tumor_mean_coverage, tumor_median_coverage, tumor_snp_overlap,
     tumor_concordance,tumor_pct_target_bases_250x, tumor_exome_coverage,
     normal_mapped_reads, normal_dedup_reads, normal_pct_target_ge_50x,
     normal_mean_coverage, normal_median_coverage, normal_snp_overlap,
     normal_concordance,normal_pct_target_bases_250x, normal_exome_coverage)
    VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
    """
    
    def to_db_value(value, numeric=False):
        if value == 'NA' or value == '' or value is None:
            return None if numeric else 'NA'
        return value

    try:
        records_inserted = 0

        with connection.cursor() as cursor:
            for row in parsed_data:
                values = (
                    to_db_value(row['run_id']),
                    to_db_value(row['case_id']),
                    to_db_value(row['tumor_id']),
                    to_db_value(row['normal_id']),
                    to_db_value(row['status']),
                    to_db_value(row['tumor_status']),
                    to_db_value(row['normal_status']),
                    int(row['tumor_mapped_reads'])          if to_db_value(row['tumor_mapped_reads'], True)          else None,
                    int(row['tumor_dedup_reads'])           if to_db_value(row['tumor_dedup_reads'], True)           else None,
                    float(row['tumor_pct_target_ge_50x'])   if to_db_value(row['tumor_pct_target_ge_50x'], True)     else None,
                    float(row['tumor_mean_coverage'])       if to_db_value(row['tumor_mean_coverage'], True)         else None,
                    float(row['tumor_median_coverage'])     if to_db_value(row['tumor_median_coverage'], True)       else None,
                    float(row['tumor_snp_overlap'])         if to_db_value(row['tumor_snp_overlap'], True)           else None,
                    to_db_value(row['tumor_concordance']),
                    float(row['tumor_pct_target_bases_250x']) if to_db_value(row['tumor_pct_target_bases_250x'], True) else None,
                    to_db_value(row['tumor_exome_coverage']),
                    int(row['normal_mapped_reads'])         if to_db_value(row['normal_mapped_reads'], True)         else None,
                    int(row['normal_dedup_reads'])          if to_db_value(row['normal_dedup_reads'], True)          else None,
                    float(row['normal_pct_target_ge_50x'])  if to_db_value(row['normal_pct_target_ge_50x'], True)    else None,
                    float(row['normal_mean_coverage'])      if to_db_value(row['normal_mean_coverage'], True)        else None,
                    float(row['normal_median_coverage'])    if to_db_value(row['normal_median_coverage'], True)      else None,
                    float(row['normal_snp_overlap'])        if to_db_value(row['normal_snp_overlap'], True)          else None,
                    to_db_value(row['normal_concordance']),
                    float(row['normal_pct_target_bases_250x']) if to_db_value(row['normal_pct_target_bases_250x'], True) else None,
                    to_db_value(row['normal_exome_coverage']),
                )
                cursor.execute(insert_query, values)
                records_inserted += 1

        connection.commit()
        print(f"✓ Successfully inserted {records_inserted} records")

    except pymysql.MySQLError as e:
        connection.rollback()
        print(f"✗ Error inserting data: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        connection.rollback()
        print(f"✗ Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


def get_table_stats(connection, table_name):
    """Get statistics about the table"""
    try:
        with connection.cursor() as cursor:
            cursor.execute(f"SELECT COUNT(*) as total FROM {table_name}")
            total = cursor.fetchone()['total']

            cursor.execute(f"""
                SELECT status, COUNT(*) as count 
                FROM {table_name} 
                GROUP BY status
            """)
            status_counts = cursor.fetchall()

            cursor.execute(f"""
                SELECT run_id, COUNT(*) as count 
                FROM {table_name} 
                GROUP BY run_id
                ORDER BY MAX(created_at) DESC
                LIMIT 5
            """)
            run_counts = cursor.fetchall()

        return {
            'total': total,
            'by_status': status_counts,
            'by_run': run_counts
        }

    except pymysql.MySQLError as e:
        print(f"⚠ Could not get table stats: {e}")
        return None


# =====================================================
# MAIN FUNCTION
# =====================================================

def main():
    """Main execution function"""
    parser = argparse.ArgumentParser(
        description='Upload paired tumor/normal QC data to existing MySQL table',
        epilog="""
Modes:
  append (default)  - Add all records to existing data
  update            - Replace data for case_ids in CSV, keep other cases

Examples:
  python %(prog)s combined_qc.csv
  python %(prog)s combined_qc.csv --update
  python %(prog)s combined_qc.csv --table other_table
        """)
    
    parser.add_argument('input_csv',
                       help='Path to input combined QC CSV file')
    parser.add_argument('--update', '-u',
                       action='store_true',
                       help='Update mode: delete existing data for case_ids in CSV before inserting')
    parser.add_argument('--table', '-t',
                       help='Table name (overrides config)')
    parser.add_argument('--yes', '-y',
                       action='store_true',
                       help='Skip confirmation prompts')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input_csv):
        print(f"✗ Error: Input file not found: {args.input_csv}", file=sys.stderr)
        sys.exit(1)
    
    table_name = args.table if args.table else TABLE_NAME
    mode = 'update' if args.update else 'append'
    
    print("=" * 80)
    print("Combined QC Upload (Tumor/Normal Paired Data)")
    print("=" * 80)
    print(f"Input File: {args.input_csv}")
    print(f"Database:   {DB_CONFIG['database']}")
    print(f"Table:      {table_name}")
    print(f"Mode:       {mode}")
    print("=" * 80)
    
    # Step 1: Parse CSV
    print("\n[Step 1/3] Parsing combined QC CSV file...")
    parsed_data = parse_combined_qc_csv(args.input_csv)
    
    if not parsed_data:
        print("✗ No data to upload")
        sys.exit(1)
    
    save_parsed_csv(parsed_data, args.input_csv)
    
    # Step 2: Connect and verify
    print("\n[Step 2/3] Connecting to database...")
    connection = create_connection(DB_CONFIG)
    
    try:
        verify_table_exists(connection, table_name)
        
        print("\nChecking for existing data...")
        duplicates = check_for_duplicates(connection, table_name, parsed_data)
        
        if duplicates:
            print(f"⚠ Found existing data for {len(duplicates)} case_id(s):")
            for case_id, count in duplicates[:5]:
                print(f"  - {case_id}: {count} records")
            if len(duplicates) > 5:
                print(f"  ... and {len(duplicates) - 5} more")
        else:
            print("✓ No duplicate case_ids found")
        
        # Step 3: Upload
        print(f"\n[Step 3/3] Uploading data...")
        
        if args.update and duplicates:
            case_ids_to_delete = [dup[0] for dup in duplicates]
            
            if not args.yes:
                print(f"\n⚠ Update mode will delete existing data for {len(case_ids_to_delete)} case_id(s)")
                response = input("Continue? Type 'yes' to confirm: ")
                if response.lower() != 'yes':
                    print("Upload cancelled.")
                    sys.exit(0)
            
            delete_existing_cases(connection, table_name, case_ids_to_delete)
        
        elif args.update and not duplicates:
            print("No existing case_ids to update, will append data")
        
        insert_parsed_data(connection, parsed_data, table_name)
        
        print("\n" + "=" * 80)
        print("TABLE STATISTICS")
        print("=" * 80)
        stats = get_table_stats(connection, table_name)
        if stats:
            print(f"\nTotal Records: {stats['total']}")
            
            if stats['by_status']:
                print("\nRecords by Status:")
                for item in stats['by_status']:
                    print(f"  {item['status']}: {item['count']}")
            
            if stats['by_run']:
                print(f"\nRecent Run IDs (top 5):")
                for item in stats['by_run']:
                    print(f"  {item['run_id']}: {item['count']} cases")
        
        print("\n" + "=" * 80)
        print("✓ Upload completed successfully!")
        print("=" * 80)
        
    finally:
        connection.close()
        print("\n✓ MySQL connection closed")


if __name__ == '__main__':
    main()