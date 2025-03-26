SET enable_progress_bar=false;
SET memory_limit='1GB';

COPY (
    SELECT 
        h.hostname || '-' || 
        (p.value->>'username') || '-' || 
        (p.value->>'command') || '-' || 
        (p.value->>'pid') AS job,
        (g.value->>'uuid')::VARCHAR AS gpuid,
        h.hostname,
        (g.value->>'name')::VARCHAR AS gputype,
        (g.value->'utilization'->>'gpu')::INTEGER AS gpu_util,
        (g.value->'memory'->>'used')::INTEGER AS gpu_mem_total_used,
        p.value->>'username' AS username,
        p.value->>'command' AS command,
        (p.value->>'pid')::INTEGER AS pid,
        (p.value->>'gpu_memory_usage')::INTEGER AS gpu_mem_usage,
        h.timestamp
    FROM read_parquet('gpu_stats_merged.parquet') h
    CROSS JOIN json_array_elements(h.json_data->'gpus') g
    CROSS JOIN json_array_elements(g.value->'processes') p
) TO 'gpu_jobs.parquet' (FORMAT PARQUET);
