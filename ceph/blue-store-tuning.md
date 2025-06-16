# Bluestore tuning

* How can i tune the blue store performance best if i have 6x24TB HDD and 8TB SSD

Here's how to optimally tune BlueStore performance for your 6x24TB HDD + 8TB SSD setup:

1. Optimal Partition Sizing
With 8TB SSD for 6 OSDs, you have plenty of space for generous allocations:

```
# Set larger DB partitions (much more than default 1%)
ceph config set osd bluestore_block_db_size 1073741824000  # 1TB per OSD

# Set larger WAL partitions  
ceph config set osd bluestore_block_wal_size 10737418240   # 10GB per OSD

# Total allocation: 6 OSDs × (1TB + 10GB) = ~6TB used, 2TB spare
```

2. BlueStore Cache Tuning

```
# Increase cache sizes for HDD OSDs
ceph config set osd bluestore_cache_size_hdd 8589934592    # 8GB cache per OSD
ceph config set osd bluestore_cache_kv_ratio 0.2           # 20% for RocksDB
ceph config set osd bluestore_cache_meta_ratio 0.8         # 80% for metadata

# Buffer sizes for large sequential writes
ceph config set osd bluestore_max_blob_size 1048576        # 1MB blobs
ceph config set osd bluestore_prefer_deferred_size 65536   # 64KB deferred writes
```


3. RocksDB Optimization

```
# Optimize RocksDB for your large DB partitions
ceph config set osd bluestore_rocksdb_options "compression=kLZ4Compression,max_write_buffer_number=32,min_write_buffer_number_to_merge=2,recycle_log_file_num=32,compaction_style=kCompactionStyleLevel,write_buffer_size=268435456,target_file_size_base=268435456,max_background_compactions=32,level0_file_num_compaction_trigger=8,level0_slowdown_writes_trigger=32,level0_stop_writes_trigger=64"
```


4. Memory and CPU Settings

```
# Increase OSD memory target (adjust based on your RAM)
ceph config set osd osd_memory_target 17179869184          # 16GB per OSD

# Optimize for your workload
ceph config set osd osd_op_num_threads_per_shard_hdd 2     # More threads for HDD
ceph config set osd osd_op_num_shards_hdd 8               # More shards
```


5. Allocation and Extent Management

```
# Optimize allocation unit size for large drives
ceph config set osd bluestore_min_alloc_size_hdd 65536     # 64KB allocation units
ceph config set osd bluestore_extent_map_shard_max_size 2048
ceph config set osd bluestore_extent_map_shard_target_size 1024
```


6. Write Path Optimization

```
# Optimize write patterns for your setup
ceph config set osd bluestore_sync_submit_transaction false
ceph config set osd bluestore_throttle_bytes 268435456     # 256MB throttle
ceph config set osd bluestore_throttle_deferred_bytes 134217728  # 128MB deferred
```

7. Deployment Command

```
# Deploy with optimized settings
ceph orch daemon add osd "$HOSTNAME:data_devices=/dev/sdb,/dev/sdc,/dev/sdd,/dev/sde,/dev/sdf,/dev/sdg,db_devices=/dev/nvme0n1,osds_per_device=1,block_db_size=1073741824000,block_wal_size=10737418240"
```

8. Post-Deployment Verification

```
# Check your partition layout
ceph-volume lvm list

# Verify configurations are applied
ceph config dump | grep bluestore

# Monitor performance
ceph osd perf
ceph osd df
```


9. Additional OS-Level Tuning

```
# Set optimal scheduler for HDDs
echo mq-deadline > /sys/block/sd*/queue/scheduler

# Set optimal scheduler for NVMe
echo none > /sys/block/nvme*/queue/scheduler

# Increase readahead for HDDs
echo 4096 > /sys/block/sd*/queue/read_ahead_kb
Key Benefits of This Configuration:
Large DB partitions (1TB each): Allows more metadata and small objects on SSD
Generous WAL (10GB each): Accommodates write bursts without spillover
Large cache: Better hit ratios for frequently accessed data
Optimized RocksDB: Better compaction and write performance
2TB spare SSD space: Available for future expansion or can be used for additional caching
```
