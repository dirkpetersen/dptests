SET enable_progress_bar=false;
SET memory_limit='1GB';

COPY (
    SELECT 
        h.hostname || '-' || 
        coalesce(p.data->>'username', '') || '-' || 
        coalesce(p.data->>'command', '') || '-' || 
        coalesce(p.data->>'pid', '') AS job,
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
    CROSS JOIN json_transform(h.json_data, '{
        "gpus": [{
            "uuid": "VARCHAR",
            "name": "VARCHAR",
            "utilization": {
                "gpu": "INTEGER"
            },
            "memory": {
                "used": "INTEGER"
            },
            "processes": [{
                "username": "VARCHAR",
                "command": "VARCHAR",
                "pid": "BIGINT",
                "gpu_memory_usage": "BIGINT"
            }]
        }]
    }') gpu_data
    CROSS JOIN unnest(gpu_data.gpus) g(data)
    CROSS JOIN unnest(g.data.processes) p(data)
) TO 'gpu_jobs.parquet' (FORMAT PARQUET);
