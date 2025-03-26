COPY (
    SELECT
        h.hostname || '-' || p.username || '-' || p.command || '-' || p.pid::VARCHAR AS job,
        g.uuid AS gpuid,
        h.hostname AS hostname,
        g.name AS gputype,
        CAST(g.utilization->>'gpu' AS INTEGER) AS gpu_util,
        CAST(g.memory->>'used' AS INTEGER) AS gpu_mem_total_used,
        p.username,
        p.command,
        CAST(p.pid AS INTEGER) AS pid,
        CAST(p.memory->>'gpu_memory_usage' AS INTEGER) AS gpu_memory_usage
    FROM 'gpu_stats_merged.parquet' AS h
    JOIN LATERAL json_transform(json(h.gpus)) AS g
        ON g.index < json_array_length(h.gpus)
    JOIN LATERAL json_transform(json(g.value->'processes')) AS p
        ON p.index < json_array_length(g.value->'processes')
) TO 'gpu_process_stats.parquet' WITH (FORMAT 'PARQUET');

