#! /usr/bin/env python3


import sqlite3
import json

# Path to the SQLite database file
db_file = "gpu_stats.sqlite"

# Name of the table in the database
table_name = "gpu_stats"

def view_table():
    try:
        # Connect to the database
        conn = sqlite3.connect(db_file)
        cursor = conn.cursor()

        # Fetch all records from the table
        cursor.execute(f"SELECT id, timestamp, hostname, json_data FROM {table_name}")
        rows = cursor.fetchall()

        # Print the records in a readable format
        for row in rows:
            print("ID:", row[0])
            print("Timestamp:", row[1])
            print("Hostname:", row[2])
            print("JSON Data:", json.dumps(json.loads(row[3]), indent=4))
            print("-" * 40)

        conn.close()
    except sqlite3.Error as e:
        print("Error while accessing the database:", e)

if __name__ == "__main__":
    view_table()

