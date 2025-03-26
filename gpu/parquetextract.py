import duckdb

# SQL query to transform nested JSON data into flattened Parquet
sql = """
COPY (
    SELECT
        hostname || '-' || username || '-' || command || '-' || pid AS job,
        json_extract(gpu.value, '$.uuid') AS gpuid,
        json_extract(json_data, '$.hostname') AS hostname,
        json_extract(gpu.value, '$.name') AS gputype,
        json_extract(gpu.value, '$.utilization.gpu') AS gpu_util,
        json_extract(gpu.value, '$.memory.used') AS gpu_mem_total_used,
        json_extract(process.value, '$.username') AS username,
        json_extract(process.value, '$.command') AS command,
        json_extract(process.value, '$.pid')::INTEGER AS pid,
        json_extract(process.value, '$.gpu_memory_usage') AS gpu_memory_usage
    FROM 'gpu_stats_merged.parquet' AS original
    CROSS JOIN LATERAL FLATTEN(input => json_extract(original.json_data, '$.gpus')) AS gpu
    CROSS JOIN LATERAL FLATTEN(input => json_extract(gpu.value, '$.processes')) AS process
) TO 'gpu_process_stats.parquet' (FORMAT PARQUET);
"""

# Execute the query with required extensions
conn = duckdb.connect()
conn.execute("INSTALL json;")  # Install JSON extension if not present
conn.execute("LOAD json;")     # Load JSON extension
conn.execute(sql)
