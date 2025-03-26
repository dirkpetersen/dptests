SET enable_progress_bar=false;
SET memory_limit='1GB';

COPY (
    SELECT 
        h.hostname || '-' || 
        COALESCE(p.username, '') || '-' || 
        COALESCE(p.command, '') || '-' || 
        COALESCE(p.pid::VARCHAR, '') AS job,
        g.uuid AS gpuid,
        h.hostname,
        g.name AS gputype,
        g.utilization->>'gpu'::INTEGER AS gpu_util,
        g.memory->>'used'::INTEGER AS gpu_mem_total_used,
        p.username,
        p.command,
        p.pid::BIGINT,
        p.gpu_memory_usage::BIGINT,
        h.timestamp
    FROM read_parquet('gpu_stats_merged.parquet') h
    CROSS JOIN json_array_elements(json_extract(h.json_data, '$.gpus')) g
    CROSS JOIN json_array_elements(json_extract(g.value, '$.processes')) p
) TO './output/gpu_jobs.parquet' (FORMAT PARQUET, COMPRESSION 'ZSTD');
