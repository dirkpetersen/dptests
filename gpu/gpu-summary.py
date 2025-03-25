import duckdb
import sys
from pathlib import Path

def main(folder_path):
    db_files = Path(folder_path).glob("*.sqlite")
    conn = duckdb.connect()
    
    # Create initial empty structure
    conn.execute("""
    CREATE TABLE combined_data AS
    SELECT 
        timestamp::TIMESTAMP AS timestamp,
        hostname,
        gpu_name,
        utilization_gpu,
        memory_used,
        username,
        command
    FROM (
        SELECT
            g.timestamp,
            g.hostname,
            j.gpus.name AS gpu_name,
            j.gpus.utilization.gpu AS utilization_gpu,
            j.gpus.memory.used AS memory_used,
            p.username,
            p.command
        FROM 
            sqlite_db.gpu_stats AS g,
            json_transform(g.json_data, '{
                "gpus": [{
                    "name": "VARCHAR",
                    "utilization": {"gpu": "INTEGER"},
                    "memory": {"used": "INTEGER"},
                    "processes": [{
                        "username": "VARCHAR", 
                        "command": "VARCHAR"
                    }]
                }]
            }') AS j,
            unnest(j.gpus) AS gpus,
            unnest(gpus.processes) AS p
    ) WHERE false
    """)

    for db_file in db_files:
        # Process each database file
        conn.execute(f"ATTACH '{str(db_file)}' AS sqlite_db (TYPE sqlite)")
        
        conn.execute("""
        INSERT INTO combined_data
        SELECT 
            timestamp::TIMESTAMP,
            hostname,
            gpu_name,
            utilization_gpu,
            memory_used,
            username,
            command
        FROM (
            SELECT
                g.timestamp,
                g.hostname,
                j.gpus.name AS gpu_name,
                j.gpus.utilization.gpu AS utilization_gpu,
                j.gpus.memory.used AS memory_used,
                p.username,
                p.command
            FROM 
                sqlite_db.gpu_stats AS g,
                json_transform(g.json_data, '{
                    "gpus": [{
                        "name": "VARCHAR",
                        "utilization": {"gpu": "INTEGER"},
                        "memory": {"used": "INTEGER"},
                        "processes": [{
                            "username": "VARCHAR", 
                            "command": "VARCHAR"
                        }]
                    }]
                }') AS j,
                unnest(j.gpus) AS gpus,
                unnest(gpus.processes) AS p
        )
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
