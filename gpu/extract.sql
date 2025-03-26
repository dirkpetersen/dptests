SET enable_progress_bar=false;
SET memory_limit='1GB';

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
        p.pid,
        CAST(p.gpu_memory_usage AS INTEGER) AS gpu_mem_usage,
        h.timestamp
    FROM read_parquet('gpu_stats_merged.parquet') h
    CROSS JOIN json_array_extract(h.json_data, '$.gpus') g
    CROSS JOIN json_array_extract(g.value, '$.processes') p
) TO 'gpu_jobs.parquet' (FORMAT PARQUET);
