SET enable_progress_bar=false;
SET memory_limit='1GB';

COPY (
    SELECT 
        h.hostname || '-' || 
        COALESCE(p.data->>'username', '') || '-' || 
        COALESCE(p.data->>'command', '') || '-' || 
        COALESCE(p.data->>'pid', '') AS job,
        g.data->>'uuid' AS gpuid,
        h.hostname,
        g.data->>'name' AS gputype,
        (g.data->'utilization'->>'gpu')::INTEGER AS gpu_util,
        (g.data->'memory'->>'used')::INTEGER AS gpu_mem_total_used,
        p.data->>'username' AS username,
        p.data->>'command' AS command,
        (p.data->>'pid')::BIGINT AS pid,
        (p.data->>'gpu_memory_usage')::BIGINT AS gpu_mem_usage,
        h.timestamp
    FROM read_parquet('gpu_stats_merged.parquet') h
    CROSS JOIN LATERAL json_array(h.json_data->'gpus') AS gu
    CROSS JOIN LATERAL json_each(gu) AS g(data)
    CROSS JOIN LATERAL json_array(g.data->'processes') AS pu
    CROSS JOIN LATERAL json_each(pu) AS p(data)
    WHERE h.json_data IS NOT NULL
) TO './output/gpu_jobs.parquet' (FORMAT PARQUET, COMPRESSION 'ZSTD');
