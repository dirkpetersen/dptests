import duckdb

# SQL query to transform nested JSON data into flattened Parquet
sql = """
COPY (
    SELECT 
        h.hostname || '-' || p.username || '-' || p.command || '-' || p.pid::VARCHAR AS job,
        g.uuid AS gpuid,
        h.hostname,
        g.name AS gputype,
        g.utilization->>'gpu' AS gpu_util,
        g.memory->>'used' AS gpu_mem_total_used,
        p.username,
        p.command,
        p.pid,
        p.memory->>'gpu_memory_usage' AS gpu_memory_usage
    FROM 'gpu_stats_merged.parquet' AS original
    JOIN LATERAL json_transform(original.json_data) AS h
    JOIN LATERAL json_transform(json(value)) AS g
        ON g.index < json_array_length(original.json_data->'$.gpus')
    JOIN LATERAL json_transform(json(g.processes)) AS p
        ON p.index < json_array_length(g.process.value)
) TO 'gpu_process_stats.parquet' (FORMAT PARQUET);
"""

# Execute the query with required extensions
conn = duckdb.connect()
conn.execute("INSTALL json")
conn.execute("LOAD json")
conn.execute(sql)
