#! /usr/bin/env python3

import os
import time
import json
import sqlite3
import subprocess

# Path to the SQLite database file
db_file = "gpu_stats.sqlite"

# Name of the table in the database
table_name = "gpu_stats"

# Path to the JSON file where gpustat output is stored
json_file = "gpu-stat-x.json"

# Ensure the database and table exist
def initialize_database():
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute(f"""
    CREATE TABLE IF NOT EXISTS {table_name} (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        timestamp TEXT,
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
    INSERT INTO {table_name} (timestamp, json_data)
    VALUES (datetime('now'), ?)
    """, (json_data,))
    conn.commit()
    conn.close()

# Ensure gpustat is installed and available
def check_gpustat():
    if not shutil.which("gpustat"):
        raise FileNotFoundError("Error: gpustat is not installed. Please install it first.")

if __name__ == "__main__":
    # Initialize the database and table
    initialize_database()

    while True:
        # Run the gpustat command and save the JSON output to a file
        try:
            result = subprocess.run(["gpustat", "--show-full-cmd", "--json"], capture_output=True, text=True)
            if result.returncode != 0:
                print("Error running gpustat:", result.stderr)
                time.sleep(900)
                continue

            with open(json_file, "w") as f:
                f.write(result.stdout)

            # Insert the JSON file as a record into the SQLite table
            insert_json_to_db(json_file)
        except Exception as e:
            print("Error:", e)

        # Wait for 15 minutes before the next iteration
        time.sleep(900)

