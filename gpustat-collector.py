#! /usr/bin/env python3

import os, sys, time, json, sqlite3, subprocess, shutil

# Path to the SQLite database file
db_file = "gpu_stats.sqlite"

# Name of the table in the database
table_name = "gpu_stats"

# Path to the JSON file where gpustat output is stored
json_file = "gpu-stat-x.json"

# Get the hostname of the system
hostname = os.uname().nodename

# Ensure the database and table exist
def initialize_database():
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute(f"""
    CREATE TABLE IF NOT EXISTS {table_name} (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        timestamp TEXT,
        hostname TEXT,
        json_data TEXT
    )
    """)
    conn.commit()
    conn.close()

# Insert JSON data into the SQLite table
def insert_json_to_db(json_path):
    with open(json_path, 'r') as f:
        json_data = f.read()

    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute(f"""
    INSERT INTO {table_name} (timestamp, hostname, json_data)
    VALUES (datetime('now'), ?, ?)
    """, (hostname, json_data))
    conn.commit()
    conn.close()

# Ensure gpustat is installed and available
def check_gpustat():
    if not shutil.which("gpustat"):
        print("Error: gpustat is not installed. Please install it first: python3 -m pip install gpustat")
        sys.exit(1)

if __name__ == "__main__":
    check_gpustat()
    # Initialize the database and table
    initialize_database()

    # Get the repeat interval from command-line arguments
    interval = None
    if len(sys.argv) > 1:
        try:
            interval = int(sys.argv[1])
        except ValueError:
            print("Error: Interval must be an integer representing seconds.")
            sys.exit(1)

    first_pass = True

    while first_pass or interval is not None:
        # Run the gpustat command and save the JSON output to a file
        try:
            result = subprocess.run(
                ["gpustat", "--show-full-cmd", "--json"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True
            )
            if result.returncode != 0:
                print("Error running gpustat:", result.stderr)
                if interval is None:
                    break
                time.sleep(interval)
                continue

            with open(json_file, "w") as f:
                f.write(result.stdout)

            # Insert the JSON file as a record into the SQLite table
            insert_json_to_db(json_file)
        except Exception as e:
            print("Error:", e)

        if interval is None:
            print(f"Stats written to table '{table_name}' in '{db_file}' !")
            break

        first_pass = False
        # Wait for the specified interval before the next iteration
        time.sleep(interval)

