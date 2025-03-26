COPY (
    SELECT
        json_extract(json_data, '$.hostname') AS hostname,
        json_extract(gpu.value, '$.uuid') AS gpuid,
        json_extract(gpu.value, '$.name') AS gputype,
        json_extract(gpu.value, '$.utilization.gpu') AS gpu_util,
        json_extract(gpu.value, '$.memory.used') AS gpu_mem_total_used,
        json_extract(process.value, '$.username') AS username,
        json_extract(process.value, '$.command') AS command,
        json_extract(process.value, '$.pid') AS pid,
        json_extract(process.value, '$.gpu_memory_usage') AS gpu_memory_usage
    FROM 'gpu_stats_merged.parquet'
    CROSS JOIN unnest(json_transform(json_extract(json_data, '$.gpus'), 'ARRAY(JSON)')) AS gpu
    LEFT JOIN unnest(json_transform(json_extract(gpu.value, '$.processes'), 'ARRAY(JSON)')) AS process
) TO 'gpu_process_stats.parquet' (FORMAT PARQUET);
