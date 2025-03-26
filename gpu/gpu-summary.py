import duckdb
import sys
from pathlib import Path

def main(folder_path):
    db_files = Path(folder_path).glob("*.sqlite")
    conn = duckdb.connect()
    
    # Load JSON extension
    conn.execute("INSTALL 'json'; LOAD 'json';")
    
    # Create empty table with explicit schema
    conn.execute("""
    CREATE TABLE combined_data (
        timestamp TIMESTAMP,
        hostname VARCHAR,
        gpu_name VARCHAR,
        utilization_gpu INTEGER,
        memory_used INTEGER,
        username VARCHAR,
        command VARCHAR
    )
    """)

    for db_file in db_files:
        # Process each database file
        conn.execute(f"ATTACH '{str(db_file)}' AS sqlite_db (TYPE SQLITE)")
        
        conn.execute("""
        INSERT INTO combined_data
        SELECT
            g.timestamp::TIMESTAMP,
            g.hostname,
            json_extract(gpu.value, '$.name')::VARCHAR AS gpu_name,
            json_extract(gpu.value, '$.utilization.gpu')::INTEGER AS utilization_gpu,
            json_extract(gpu.value, '$.memory.used')::INTEGER AS memory_used,
            json_extract(process.value, '$.username')::VARCHAR AS username,
            json_extract(process.value, '$.command')::VARCHAR AS command
        FROM 
            sqlite_db.gpu_stats AS g,
            unnest(json_extract(g.json_data, '$.gpus')) AS gpu,
            unnest(json_extract(gpu.value, '$.processes')) AS process
        """)
        
        conn.execute("DETACH sqlite_db")

    # Export to Parquet
    conn.execute("COPY combined_data TO 'gpu_metrics.parquet' (FORMAT PARQUET)")
    conn.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python gpu-summary.py <folder_path>")
        sys.exit(1)
    main(sys.argv[1])
