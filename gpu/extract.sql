SET enable_progress_bar=false;
SET memory_limit='1GB';

-- Create output directory
CREATE TABLE IF NOT EXISTS duckdb_fs_directories(path VARCHAR);
INSERT OR REPLACE INTO duckdb_fs_directories SELECT './output' WHERE NOT EXISTS (SELECT 1 FROM duckdb_fs_directories WHERE path = './output');

-- Export query results
COPY (
    SELECT 
        h.hostname || '-' || 
        COALESCE(p.value->>'username', '') || '-' || 
        COALESCE(p.value->>'command', '') || '-' || 
        COALESCE(p.value->>'pid', '') AS job,
        g.value->>'uuid' AS gpuid,
        h.hostname,
        g.value->>'name' AS gputype,
        (g.value->'utilization'->>'gpu')::INTEGER AS gpu_util,
        (g.value->'memory'->>'used')::INTEGER AS gpu_mem_total_used,
        p.value->>'username' AS username,
        p.value->>'command' AS command,
        (p.value->>'pid')::BIGINT AS pid,
        (p.value->>'gpu_memory_usage')::BIGINT AS gpu_mem_usage,
        h.timestamp
    FROM read_parquet('gpu_stats_merged.parquet') h
    CROSS JOIN LATERAL json_extract(h.json_data, '$.gpus') gpus
    CROSS JOIN LATERAL json_each(gpus) g
    CROSS JOIN LATERAL json_extract(g.value, '$.processes') procs
    CROSS JOIN LATERAL json_each(procs) p
    WHERE h.json_data IS NOT NULL
) TO './output/gpu_jobs.parquet' (FORMAT PARQUET, COMPRESSION 'ZSTD');
