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
import re

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
# Safely validate and quote MySQL table names before using them in SQL
def quote_identifier(identifier):
    """
    Safely quote a MySQL identifier such as a table name.

    SQL parameters protect values, but not table names. Since table_name is
    interpolated into SQL strings, restrict it to ordinary identifier characters.
    """
    if not re.fullmatch(r"[A-Za-z0-9_]+", identifier):
        raise ValueError(f"Unsafe SQL identifier: {identifier}")

    return f"`{identifier}`"


# NEW FUNC: Build the logical duplicate keys for control records.
def get_control_upload_keys(parsed_data):
    """
    Return the logical upload keys for this controls file.

    A duplicate control row should usually mean the same sample in the same run,
    not merely any row from the same run_id.
    """
    keys = set()
    INVALID_SAMPLE_NAMES = {"", "NA", "N/A", "Not Found", "NULL", "None"}

    for row in parsed_data:
        run_id = str(row.get("run_id", "")).strip()
        sample_name = str(row.get("sample_name", "")).strip()

        if run_id in INVALID_SAMPLE_NAMES:
            continue

        if sample_name in INVALID_SAMPLE_NAMES:
            continue

        keys.add((run_id, sample_name))

    return sorted(keys)


def parse_qc_csv(input_file):
    """Parse the QC CSV file and extract required columns"""
    column_mapping = {
        "Run_ID": "run_id",
        "RunQC_Sample": "sample_name",
        "Overall_QC": "sample_status",
        "Mapped Reads QC_Value": "mapped_reads",
        "Overall Coverage_Value": "coverage",
        "Negative control QC_Value": "negative_control",
        "Positive control QC_Value": "positive_control",
    }

    parsed_data = []

    try:
        with open(input_file, "r") as infile:
            reader = csv.DictReader(infile)

            missing = [col for col in column_mapping if col not in (reader.fieldnames or [])]
            if missing:
                print(f"✗ Error: Missing expected columns: {missing}", file=sys.stderr)
                print(f"  Found columns: {reader.fieldnames}", file=sys.stderr)

            for row in reader:
                parsed_row = {}
                for old_name, new_name in column_mapping.items():
                    value = row.get(old_name, "").strip()
                    parsed_row[new_name] = "NA" if value == "" else value
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
        output_file = f"{name_without_ext}_parsed.csv"

        if parsed_data:
            fieldnames = [
                "run_id",
                "sample_name",
                "sample_status",
                "mapped_reads",
                "coverage",
                "negative_control",
                "positive_control",
            ]

            with open(output_file, "w", newline="") as outfile:
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
            host=config["host"],
            user=config["user"],
            password=config["password"],
            database=config["database"],
            port=config.get("port", 3306),
            cursorclass=pymysql.cursors.DictCursor,
        )
        print("✓ Connected to MySQL Server")
        print(f"✓ Using database: {config['database']}")
        return connection

    except pymysql.MySQLError as e:
        print(f"✗ Error connecting to MySQL: {e}", file=sys.stderr)
        sys.exit(1)


# Table check uses a parameterized value instead of directly from table_name
def verify_table_exists(connection, table_name):
    """Verify table exists - exit if it doesn't."""
    try:
        with connection.cursor() as cursor:
            cursor.execute("SHOW TABLES LIKE %s", (table_name,))
            result = cursor.fetchone()

        if not result:
            print(f"✗ Error: Table '{table_name}' does not exist", file=sys.stderr)
            print("This script only works with existing tables.")
            sys.exit(1)

        print(f"✓ Table '{table_name}' exists")

    except pymysql.MySQLError as e:
        print(f"✗ Error checking table: {e}", file=sys.stderr)
        sys.exit(1)


# Duplicate check targets specific run_id + sample_name pairs instead of just run_id
def check_for_duplicates(connection, table_name, parsed_data):
    """
    Check whether records already exist for the same run_id + sample_name pairs.

    This avoids false positives where the same run already exists but the
    uploaded control sample is new.
    """
    try:
        upload_keys = get_control_upload_keys(parsed_data)

        if not upload_keys:
            return []

        safe_table_name = quote_identifier(table_name)
        placeholders = ", ".join(["(%s, %s)"] * len(upload_keys))
        params = []

        for run_id, sample_name in upload_keys:
            params.extend([run_id, sample_name])

        query = f"""
            SELECT run_id, sample_name, COUNT(*) AS count
            FROM {safe_table_name}
            WHERE (run_id, sample_name) IN ({placeholders})
            GROUP BY run_id, sample_name
        """

        with connection.cursor() as cursor:
            cursor.execute(query, params)
            rows = cursor.fetchall()

        return [(row["run_id"], row["sample_name"], row["count"]) for row in rows]

    except Exception as e:
        print(f"⚠ Could not check for duplicates: {e}")
        return []


# Convert integer-like values safely before database insert
def to_db_int(value):
    """Convert integer-like values to int, returning None for missing values."""
    if value is None:
        return None

    value = str(value).strip()

    if value in {"", "NA", "N/A", "NULL", "None"}:
        return None

    return int(float(value.replace(",", "")))


# Convert float-like values safely before database insert
def to_db_float(value):
    """Convert float-like values to float, returning None for missing values."""
    if value is None:
        return None

    value = str(value).strip()

    if value in {"", "NA", "N/A", "NULL", "None"}:
        return None

    return float(value.replace(",", "").replace("%", ""))


# Update/delete logic targets specific run_id + sample_name pairs instead of just run_id
def delete_existing_controls(connection, table_name, record_keys):
    """
    Delete existing rows for specific run_id + sample_name pairs.

    This prevents update mode from deleting all controls from a run when only
    selected control samples are being replaced.
    """
    try:
        if not record_keys:
            print("✓ No existing control records to delete")
            return

        safe_table_name = quote_identifier(table_name)
        placeholders = ", ".join(["(%s, %s)"] * len(record_keys))
        params = []

        for run_id, sample_name in record_keys:
            params.extend([run_id, sample_name])

        query = f"""
            DELETE FROM {safe_table_name}
            WHERE (run_id, sample_name) IN ({placeholders})
        """

        with connection.cursor() as cursor:
            cursor.execute(query, params)
            deleted = cursor.rowcount

        connection.commit()
        print(
            f"✓ Deleted {deleted} existing record(s) for {len(record_keys)} run/sample pair(s)"
        )

    except pymysql.MySQLError as e:
        connection.rollback()
        print(f"✗ Error deleting data: {e}", file=sys.stderr)
        sys.exit(1)


# Insert uses a quoted table name and safer numeric conversion helpers
def insert_parsed_data(connection, parsed_data, table_name):
    """Insert parsed data into the database."""
    INVALID_SAMPLE_NAMES = {"", "NA", "N/A", "Not Found", "NULL", "None"}

    safe_table_name = quote_identifier(table_name)

    insert_query = f"""
    INSERT INTO {safe_table_name}
    (run_id, sample_name, sample_status, mapped_reads, coverage,
     negative_control, positive_control)
    VALUES (%s, %s, %s, %s, %s, %s, %s)
    """

    try:
        records_inserted = 0

        with connection.cursor() as cursor:
            for row in parsed_data:
                if row["sample_name"].strip() in INVALID_SAMPLE_NAMES:
                    print(f"⚠ Skipping row with invalid sample_name: {row}")
                    continue
                values = (
                    row["run_id"],
                    row["sample_name"],
                    row["sample_status"],
                    to_db_int(row["mapped_reads"]),
                    to_db_float(row["coverage"]),
                    row["negative_control"],
                    row["positive_control"],
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


# Table stats no longer assumes created_at exists
def get_table_stats(connection, table_name):
    """Get statistics about the table."""
    try:
        safe_table_name = quote_identifier(table_name)

        with connection.cursor() as cursor:
            cursor.execute(f"SELECT COUNT(*) AS total FROM {safe_table_name}")
            total = cursor.fetchone()["total"]

            cursor.execute(f"""
                SELECT sample_status, COUNT(*) AS count
                FROM {safe_table_name}
                GROUP BY sample_status
            """)
            status_counts = cursor.fetchall()

            cursor.execute(f"""
                SELECT run_id, COUNT(*) AS count
                FROM {safe_table_name}
                GROUP BY run_id
                ORDER BY run_id DESC
                LIMIT 5
            """)
            run_counts = cursor.fetchall()

        return {
            "total": total,
            "by_status": status_counts,
            "by_run": run_counts,
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
        description="Upload QC data to existing MySQL table (safe - no table creation/clearing)",
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
        """,
    )

    parser.add_argument("input_csv", help="Path to input QC CSV file")
    parser.add_argument("--update", "-u",
        action="store_true",
        help="Update mode: delete existing data for run_ids in CSV before inserting",
    )
    parser.add_argument("--table", "-t", help="Table name (overrides config)")
    parser.add_argument("--yes", "-y", action="store_true",
                        help="Skip confirmation prompts"
    )

    args = parser.parse_args()

    if not os.path.exists(args.input_csv):
        print(f"✗ Error: Input file not found: {args.input_csv}", file=sys.stderr)
        sys.exit(1)

    table_name = args.table if args.table else CONTROL_TABLE_NAME
    mode = "update" if args.update else "append"

    print("=" * 80)
    print("Safe QC Upload to Existing Table")
    print("=" * 80)
    print(f"Input File: {args.input_csv}")
    print(f"Database:   {DB_CONFIG['database']}")
    print(f"Table:      {table_name}")
    print(f"Mode:       {mode}")
    print("=" * 80)

    # Step 1: Parse the CSV ----------------------------------------------
    print("\n[Step 1/3] Parsing QC CSV file...")
    parsed_data = parse_qc_csv(args.input_csv)

    if not parsed_data:
        print("✗ No data to upload")
        sys.exit(1)

    save_parsed_csv(parsed_data, args.input_csv)

    # Step 2: Connect and verify ------------------------------------------
    print("\n[Step 2/3] Connecting to database...")
    connection = create_connection(DB_CONFIG)

    try:
        verify_table_exists(connection, table_name)

        # Duplicate display reports exact run_id + sample_name pairs
        print("\nChecking for existing data...")
        duplicates = check_for_duplicates(connection, table_name, parsed_data)

        if duplicates:
            print(f"⚠ Found existing data for {len(duplicates)} run/sample pair(s):")
            for run_id, sample_name, count in duplicates[:5]:
                print(
                    f"  - run_id={run_id}, sample_name={sample_name}: {count} record(s)"
                )

            if len(duplicates) > 5:
                print(f"  ... and {len(duplicates) - 5} more")
        else:
            print("✓ No duplicate run/sample pairs found")

        # Step 3: Upload ----------------------------------------------
        print("\n[Step 3/3] Uploading data...")

        # Update/delete logic targets specific run_id + sample_name pairs 
        # instead of just run_id.
        if args.update and duplicates:
            record_keys_to_delete = [
                (run_id, sample_name) for run_id, sample_name, count in duplicates
            ]

            if not args.yes:
                print(
                    f"\n⚠ Update mode will delete existing data for "
                    f"{len(record_keys_to_delete)} run/sample pair(s)"
                )
                response = input("Continue? Type 'yes' to confirm: ")

                if response.lower() != "yes":
                    print("Upload cancelled.")
                    sys.exit(0)

            delete_existing_controls(connection, table_name, record_keys_to_delete)

        elif args.update and not duplicates:
            print("No existing run/sample pairs to update, will append data")

        insert_parsed_data(connection, parsed_data, table_name)

        print("\n" + "=" * 80)
        print("TABLE STATISTICS")
        print("=" * 80)
        stats = get_table_stats(connection, table_name)
        if stats:
            print(f"\nTotal Records: {stats['total']}")

            if stats["by_status"]:
                print("\nRecords by Status:")
                for item in stats["by_status"]:
                    print(f"  {item['sample_status']}: {item['count']}")

            if stats["by_run"]:
                print("\nRecent Run IDs (top 5):")
                for item in stats["by_run"]:
                    print(f"  {item['run_id']}: {item['count']} samples")

        print("\n" + "=" * 80)
        print("✓ Upload completed successfully!")
        print("=" * 80)

    finally:
        connection.close()
        print("\n✓ MySQL connection closed")


if __name__ == "__main__":
    main()
