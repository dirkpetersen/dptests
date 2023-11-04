#! /usr/bin/env python3

import pandas, sqlite3

# File paths
csv_file_path = './db/CEDAR.csv'
sqlite_db_path = './db/CEDAR.db'

# Read CSV file into dataframe and copy the dataframe to SQLite
df = pandas.read_csv(csv_file_path)

# Connect to SQLite database (it will be created if it doesn’t exist)
conn = sqlite3.connect(sqlite_db_path)

# Insert data
df.to_sql('meta_table', conn, index=False, if_exists='replace')

# Commit and close
conn.commit()
conn.close()


