#!/bin/bash

# Path to the DuckDB database file
db_file="gpu_stats.duckdb"

# Name of the table in the database
table_name="gpu_stats"

# Path to the JSON file where gpustat output is stored
json_file="gpu-stat-x.json"

# Ensure DuckDB is installed and available
if ! command -v duckdb &> /dev/null; then
    echo "Error: DuckDB is not installed. Please install it first."
    exit 1
fi

# Ensure gpustat is installed and available
if ! command -v gpustat &> /dev/null; then
    echo "Error: gpustat is not installed. Please install it first."
    exit 1
fi

while true; do
    # Run the gpustat command and save the JSON output to a file
    gpustat --show-full-cmd --json > "$json_file"

    # Insert the JSON file as a record into the DuckDB table
    duckdb "$db_file" <<SQL
    INSERT INTO $table_name
    SELECT *
    FROM read_json_auto('$json_file');
SQL

    # Wait for 15 minutes before the next iteration
    sleep 900
done

