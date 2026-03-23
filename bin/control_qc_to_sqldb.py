#!/usr/bin/env python3
"""
Safe QC Upload - For existing database and table
Only appends or updates data, never creates or clears tables
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
    from db_config import DB_CONFIG, CONTROL_TABLE_NAME
    print("✓ Using configuration from db_config.py")
except ImportError:
    print("✗ Error: db_config.py not found")
    print("Please create db_config.py from db_config_template.py")
    sys.exit(1)


# =====================================================
# PARSING FUNCTIONS
# =====================================================

def parse_qc_csv(input_file):
    """Parse the QC CSV file and extract required columns"""
    column_mapping = {
        'Run_ID': 'run_id',
        'RunQC_Sample': 'sample_name',
        'Overall_QC': 'sample_status',
        'Mapped Reads QC_Value': 'mapped_reads',
        'Overall Coverage_Value': 'coverage',
        'Negative control QC_Value': 'negative_control',
        'Positive control QC_Value': 'positive_control'
    }
    
    parsed_data = []
    
    try:
        with open(input_file, 'r') as infile:
            reader = csv.DictReader(infile)
            
            for row in reader:
                parsed_row = {}
                for old_name, new_name in column_mapping.items():
                    value = row.get(old_name, '').strip()
                    # Convert empty values to 'NA'
                    if value == '' or value is None:
                        parsed_row[new_name] = 'NA'
                    else:
                        parsed_row[new_name] = value
                
                parsed_data.append(parsed_row)
        
        print(f"✓ Parsed {len(parsed_data)} records from {input_file}")
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
            fieldnames = ['run_id', 'sample_name', 'sample_status', 'mapped_reads',
                         'coverage', 'negative_control', 'positive_control']
            
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
        run_ids = set(row['run_id'] for row in parsed_data if row['run_id'])
        duplicates = []

        with connection.cursor() as cursor:
            for run_id in run_ids:
                cursor.execute(
                    f"SELECT COUNT(*) as count FROM {table_name} WHERE run_id = %s",
                    (run_id,)
                )
                count = cursor.fetchone()['count']
                if count > 0:
                    duplicates.append((run_id, count))

        return duplicates

    except Exception as e:
        print(f"⚠ Could not check for duplicates: {e}")
        return []


def delete_existing_runs(connection, table_name, run_ids):
    """Delete existing data for specific run_ids"""
    try:
        placeholders = ','.join(['%s'] * len(run_ids))
        query = f"DELETE FROM {table_name} WHERE run_id IN ({placeholders})"

        with connection.cursor() as cursor:
            cursor.execute(query, list(run_ids))
            deleted = cursor.rowcount

        connection.commit()
        print(f"✓ Deleted {deleted} existing records for {len(run_ids)} run_id(s)")

    except pymysql.MySQLError as e:
        connection.rollback()
        print(f"✗ Error deleting data: {e}", file=sys.stderr)
        sys.exit(1)


def insert_parsed_data(connection, parsed_data, table_name):
    """Insert parsed data into the database"""
    insert_query = f"""
    INSERT INTO {table_name} 
    (run_id, sample_name, sample_status, mapped_reads, coverage,
     negative_control, positive_control)
    VALUES (%s, %s, %s, %s, %s, %s, %s)
    """
    
    try:
        records_inserted = 0

        with connection.cursor() as cursor:
            for row in parsed_data:
                values = (
                    row['run_id'],
                    row['sample_name'],
                    row['sample_status'],
                    int(row['mapped_reads']) if row['mapped_reads'] != 'NA' else None,
                    float(row['coverage']) if row['coverage'] != 'NA' else None,
                    row['negative_control'],
                    row['positive_control']
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
                SELECT sample_status, COUNT(*) as count 
                FROM {table_name} 
                GROUP BY sample_status
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
        description='Upload QC data to existing MySQL table (safe - no table creation/clearing)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Modes:
  append (default)  - Add all records to existing data
  update            - Replace data for run_ids in CSV, keep other runs

Examples:
  # Append new data
  python %(prog)s controls_qc.csv

  # Update existing runs
  python %(prog)s controls_qc.csv --update

  # Use different table
  python %(prog)s data.csv --table other_table
        """)
    
    parser.add_argument('input_csv',
                       help='Path to input QC CSV file')
    parser.add_argument('--update', '-u',
                       action='store_true',
                       help='Update mode: delete existing data for run_ids in CSV before inserting')
    parser.add_argument('--table', '-t',
                       help='Table name (overrides config)')
    parser.add_argument('--yes', '-y',
                       action='store_true',
                       help='Skip confirmation prompts')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input_csv):
        print(f"✗ Error: Input file not found: {args.input_csv}", file=sys.stderr)
        sys.exit(1)
    
    table_name = args.table if args.table else CONTROL_TABLE_NAME
    mode = 'update' if args.update else 'append'
    
    print("=" * 80)
    print("Safe QC Upload to Existing Table")
    print("=" * 80)
    print(f"Input File: {args.input_csv}")
    print(f"Database:   {DB_CONFIG['database']}")
    print(f"Table:      {table_name}")
    print(f"Mode:       {mode}")
    print("=" * 80)
    
    # Step 1: Parse the CSV
    print("\n[Step 1/3] Parsing QC CSV file...")
    parsed_data = parse_qc_csv(args.input_csv)
    
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
            print(f"⚠ Found existing data for {len(duplicates)} run_id(s):")
            for run_id, count in duplicates[:5]:
                print(f"  - {run_id}: {count} records")
            if len(duplicates) > 5:
                print(f"  ... and {len(duplicates) - 5} more")
        else:
            print("✓ No duplicate run_ids found")
        
        # Step 3: Upload
        print(f"\n[Step 3/3] Uploading data...")
        
        if args.update and duplicates:
            run_ids_to_delete = [dup[0] for dup in duplicates]
            
            if not args.yes:
                print(f"\n⚠ Update mode will delete existing data for {len(run_ids_to_delete)} run_id(s)")
                response = input("Continue? Type 'yes' to confirm: ")
                if response.lower() != 'yes':
                    print("Upload cancelled.")
                    sys.exit(0)
            
            delete_existing_runs(connection, table_name, run_ids_to_delete)
        
        elif args.update and not duplicates:
            print("No existing run_ids to update, will append data")
        
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
                    print(f"  {item['sample_status']}: {item['count']}")
            
            if stats['by_run']:
                print(f"\nRecent Run IDs (top 5):")
                for item in stats['by_run']:
                    print(f"  {item['run_id']}: {item['count']} samples")
        
        print("\n" + "=" * 80)
        print("✓ Upload completed successfully!")
        print("=" * 80)
        
    finally:
        connection.close()
        print("\n✓ MySQL connection closed")


if __name__ == '__main__':
    main()