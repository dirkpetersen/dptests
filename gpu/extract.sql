COPY (
    SELECT
        json_extract(main.json_data, '$.hostname') AS hostname,
        json_extract(gpu.value, '$.uuid') AS gpuid,
        json_extract(gpu.value, '$.name') AS gputype,
        json_extract(gpu.value, '$.utilization.gpu')::INTEGER AS gpu_util,
        json_extract(gpu.value, '$.memory.used')::INTEGER AS gpu_mem_total_used,
        json_extract(process.value, '$.username') AS username,
        json_extract(process.value, '$.command') AS command,
        json_extract(process.value, '$.pid')::INTEGER AS pid,
        json_extract(process.value, '$.gpu_memory_usage')::INTEGER AS gpu_memory_usage
    FROM 'gpu_stats_merged.parquet' AS main
    CROSS JOIN LATERAL unnest(json_transform(json_extract(main.json_data, '$.gpus'), 'JSON[]')) AS gpu
    LEFT JOIN LATERAL unnest(json_transform(json_extract(gpu.value, '$.processes'), 'JSON[]')) AS process
) AS processed_data TO 'gpu_process_stats.parquet' WITH (FORMAT 'PARQUET');
