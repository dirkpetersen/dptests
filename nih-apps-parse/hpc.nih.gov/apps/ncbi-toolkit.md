

document.querySelector('title').textContent = 'ncbi-toolkit on Biowulf';
ncbi-toolkit on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
[List of executables](#list)
 |


The [NCBI C++ Toolkit](http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/) is a set of executables and
libraries for a multitude of sequence analysis functions.


These executables have been compiled and made available.


Documentation
* [The NCBI C++ Toolkit Book](https://ncbi.github.io/cxx-toolkit/)


Many of the executables have help functions. These can be displayed with the **-help** option:



```

**$ fastq-dump -help**

Usage:
  fastq-dump [options] 
 fastq-dump [options] [ -A ] 

INPUT
 -A|--accession  Replaces accession derived from  in 
 filename(s) and deflines (only for single 
 table dump) 
 --table  Table name within cSRA object, default is 
 "SEQUENCE" 

PROCESSING

Read Splitting Sequence data may be used in raw form or
 split into individual reads
 --split-spot Split spots into individual reads 

Full Spot Filters Applied to the full spot independently
 of --split-spot
 -N|--minSpotId  Minimum spot id 
 -X|--maxSpotId  Maximum spot id 
 --spot-groups <[list]> Filter by SPOT\_GROUP (member): name[,...] 
 -W|--clip Apply left and right clips 

... etc ...

```

Important Notes
* Module Name: ncbi-toolkit (see [the modules page](/apps/modules.html) for more information)



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load ncbi-toolkit**
[user@cn3144 ~]$ **gi2taxid -gi 36209385**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. ncbi-toolkit.sh). For example:



```
#!/bin/bash
module load ncbi-toolkit

# NOTE: This is merely a test to see that ncbi-toolkit runs correctly.
# This example may not be rational or useful.

# Create a nucleotide blast database suitable for the ncbi-toolkit version.
# In this example, we extract the top 1,000,000 lines from nt.fas.

head -1000000 /fdb/fastadb/nt.fas > nt_1M.fas
makeblastdb -in nt_1M.fas -dbtype nucl

# Now run Repeat Masker blast against this database.

rmblastn -query gi_255958152.nt.fas -db nt_1M.fas -gapopen 3 -gapextend 3
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] ncbi-toolkit.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. ncbi-toolkit.swarm). For example:



```

igblastn -db mydb -query seq1.fas -out seq1.out
igblastn -db mydb -query seq2.fas -out seq2.out
igblastn -db mydb -query seq3.fas -out seq3.out
igblastn -db mydb -query seq4.fas -out seq4.out

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f ncbi-toolkit.swarm [-g #] [-t #] --module ncbi-toolkit
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module ncbi-toolkit Loads the ncbi-toolkit module for each subjob in the swarm 
 | |
 | |
 | |


List of executables (as in v21.0.0)
-----------------------------------



```

aalookup_unit_test           lmdb_test4                      test_grid_worker
aascan_unit_test             lmdb_test5                      test_hash
ace2asn                      lmdbxx_sample                   test_hgvs_parser
adapter_search               localfinder                     test_hgvs_reader
agpconvert                   logs_replay                     test_histogram_binning
agp_count                    logs_splitter                   test_html
agp_renumber                 magicblast                      test_ic_client
agp_validate                 magicblast_unit_test            test_id1_client
agp_val_test                 makeblastdb                     test_id_mux
align_filter_unit_test       makembindex                     test_image
align_format_unit_test       makeprofiledb                   test_interprocess_lock
aln_build                    mapper_unit_test                test_json_over_uttp
alnmgr_sample                merge_tree                      test_jsonwrapp
alnmrg                       mirror_test                     test_jsonwrapp10
aln_test                     msa2pssm_unit_test              test_lds2
alnvwr                       multi_command                   test_limited_map
annotwriter                  multireader                     test_line_reader
asn2asn                      mysql_lang                      test_logrotate
asn2fasta                    ncbi_applog                     test_math
asn2flat                     ncbi_dblb_cli                   test_message_mt
asn_assign                   ncbi_encrypt                    test_metaphone
asn_cache_test               ncfetch.cgi                     test_multipart.cgi
asn_cleanup                  netcache_cgi_sample.cgi         test_ncbiargs
asndisc                      netcache_client_sample          test_ncbiargs_sample
asniotest                    netcache_control                test_ncbi_buffer
asn_sample                   netcached                       test_ncbicfg
asnvalidate                  netschedule_client_sample       test_ncbi_clog_mt
asnwalk_read                 netschedule_control             test_ncbi_clog_mt_ctx
asnwalk_type                 netscheduled                    test_ncbi_clog_templates
asnwalk_write                netschedule_node_sample         test_ncbicntr
autodef_demo                 netstoraged                     test_ncbi_config
bam2graph                    netstorage_gc                   test_ncbi_conn
bamgraph_test                ngalign_app                     test_ncbi_conn_stream
bamindex_test                nmer_repeats                    test_ncbi_conn_stream_mt
bam_test                     ns_loader                       test_ncbi_connutil_hit
basic_sample                 ns_remote_job_control           test_ncbi_connutil_misc
basic_sample_lib_test        ns_submit_remote_job            test_ncbi_core
bdb_demo1                    ntlookup_unit_test              test_ncbidiag_f_mt
bdb_demo2                    ntscan_unit_test                test_ncbidiag_mt
bdb_demo3                    nw_aligner                      test_ncbidiag_p
bdb_demo4                    objects_sample                  test_ncbi_disp
bdb_dump                     objextract                      test_ncbidll
bdb_env_keeper               objmgr_demo                     test_ncbi_download
bdbloader_unit_test          objmgr_sample                   test_ncbi_dsock
bdbtest                      odbc95_array                    test_ncbiexec
bdbtest_split                odbc95_array_error              test_ncbiexpt
biosample_chk                odbc95_array_out                test_ncbi_fast
bioseq_edit_sample           odbc95_attributes               test_ncbifile
bl2seq_unit_test             odbc95_binary_test              test_ncbi_file_connector
blast_dataloader_unit_test   odbc95_blob1                    test_ncbi_ftp_connector
blastdb_aliastool            odbc95_cancel                   test_ncbi_ftp_download
blastdbcheck                 odbc95_closestmt                test_ncbi_heapmgr
blastdbcmd                   odbc95_compute                  test_ncbi_hmac
blastdb_convert              odbc95_connect                  test_ncbi_http_connector
blastdbcp                    odbc95_connect2                 test_ncbi_http_get
blastdb_format_unit_test     odbc95_const_params             test_ncbi_http_upload
blast_demo                   odbc95_convert_error            test_ncbi_iprange
blastdiag_unit_test          odbc95_copydesc                 test_ncbi_ipv6
blastengine_unit_test        odbc95_cursor1                  test_ncbi_limits
blastextend_unit_test        odbc95_cursor2                  test_ncbi_linkerd
blastfilter_unit_test        odbc95_cursor3                  test_ncbi_linkerd_mt
blast_formatter              odbc95_cursor4                  test_ncbi_memory_connector
blast_format_unit_test       odbc95_cursor5                  test_ncbimime
blasthits_unit_test          odbc95_cursor6                  test_ncbi_namedpipe
blastinput_demo              odbc95_cursor7                  test_ncbi_namedpipe_connector
blastinput_unit_test         odbc95_data                     test_ncbi_namerd
blastn                       odbc95_date                     test_ncbi_namerd_mt
blastoptions_unit_test       odbc95_descrec                  test_ncbi_null
blastp                       odbc95_describecol              test_ncbi_os_unix
blast_sample                 odbc95_describecol2             test_ncbi_pipe
blast_services_unit_test     odbc95_earlybind                test_ncbi_pipe_connector
blastsetup_unit_test         odbc95_error                    test_ncbi_process
blastsrainput_unit_test      odbc95_freeclose                test_ncbi_rate_monitor
blast_tabular_unit_test      odbc95_funccall                 test_ncbireg_mt
blast_unit_test              odbc95_genparams                test_ncbi_rwstream
blastx                       odbc95_getdata                  test_ncbi_sendmail
blobreader                   odbc95_hidden                   test_ncbi_server_info
blobrwd                      odbc95_insert_speed             test_ncbi_service
blobrws                      odbc95_lang_error               test_ncbi_service_connector
blobwriter                   odbc95_long_error               test_ncbi_service_cxx_mt
bma_refiner                  odbc95_mars1                    test_ncbi_socket
bm_sparse_sample             odbc95_moreandcount             test_ncbi_socket_connector
bss_info                     odbc95_norowset                 test_ncbistr
cache_demo                   odbc95_paramcore                test_ncbi_system
cache_index_copy             odbc95_params                   test_ncbi_table
cddalignview                 odbc95_peter                    test_ncbithr
cgi2rcgi                     odbc95_prepare_results          test_ncbithr_native
cgi_io_test                  odbc95_prepare_warn             test_ncbitime
cgi_redirect                 odbc95_prepclose                test_ncbitime_mt
cgi_sample.cgi               odbc95_preperror                test_ncbi_tree
cgi_session_sample.cgi       odbc95_print                    test_ncbi_trigger
cgitest                      odbc95_putdata                  test_ncbi_url
cgi_tunnel2grid.cgi          odbc95_raiserror                test_ncbiutil
cidtool                      odbc95_rebindpar                test_nc_stress
clusterer                    odbc95_rowset                   test_nc_stress_pubmed
cobalt                       odbc95_rpc                      test_netcache_api
cobalt_unit_test             odbc95_scroll                   test_netschedule_client
collection_scores_unit_test  odbc95_stats                    test_netschedule_crash
compart                      odbc95_t0001                    test_netschedule_node
compartp                     odbc95_t0002                    test_netschedule_stress
concat_seqentries            odbc95_t0003                    test_netservice_params
convert2blastmask            odbc95_t0004                    test_netstorage
convert_seq                  odbc95_tables                   test_nsstorage
conv_image                   odbc95_test64                   test_objmgr
coretest                     odbc95_testodbc                 test_objmgr_basic
cpgdemo                      odbc95_timeout                  test_objmgr_data
csra_test_mt                 odbc95_timeout2                 test_objmgr_data_mt
ct95_array_bind              odbc95_timeout3                 test_objmgr_gbloader
ct95_blk_in                  odbc95_timeout4                 test_objmgr_gbloader_mt
ct95_blk_in2                 odbc95_transaction              test_objmgr_mem
ct95_blk_out                 odbc95_transaction2             test_objmgr_mt
ct95_cancel                  odbc95_type                     test_objmgr_sv
ct95_connect_fail            odbc95_typeinfo                 test_objstore
ct95_cs_config               odbc95_utf8                     test-odbc95
ct95_cs_diag                 odbc95_utf8_2                   test_param_mt
ct95_ct_cursor               odbc95_warning                  test_plugins
ct95_ct_cursors              odbc95_wchar                    test_porter_stemming
ct95_ct_diagall              omssa2pepXML                    test_queue_mt
ct95_ct_diagclient           omssacl                         test_random
ct95_ct_diagserver           omssamerge                      test_range_coll
ct95_ct_dynamic              optionshandle_unit_test         test_rangemap
ct95_ct_options              pacc                            test_reader_gicache
ct95_datafmt                 phiblast_unit_test              test_reader_id1
ct95_get_send_data           phytree_calc_unit_test          test-reference_allele_fix
ct95_lang_ct_param           phytree_format_unit_test        test_regexp
ct95_rpc_ct_param            prelimsearch_unit_test          test_relloc
ct95_rpc_ct_setparam         prime_cache                     test_request_control
ct95_rpc_fail                project_tree_builder            test_resize_iter
ct95_t0001                   proteinkmer_unit_test           test_resource_info
ct95_t0002                   prot_match                      test_row_reader
ct95_t0003                   psibl2seq_unit_test             test_row_reader_excel_csv
ct95_t0004                   psiblast                        test_row_reader_iana_csv
ct95_t0005                   psiblast_iteration_unit_test    test_row_reader_iana_tsv
ct95_t0006                   psiblast_unit_test              test_row_reader_ncbi_tsv
ct95_t0007                   pssmcreate_unit_test            test_row_reader_performance
ct95_t0008                   pssmenginefreqratios_unit_test  test_scheduler
ct95_t0009                   pub_report                      test_scoremat
ctl_lang_ftds64              python_ncbi_dbapi_test          test_semaphore_mt
ctl_lang_ftds95              querydata_unit_test             test_seq_entry_ci
ctl_sp_databases_ftds        queryinfo_unit_test             test_seqio
ctl_sp_who_ftds64            rcgi_sample                     test_seqmap_switch
ctl_sp_who_ftds95            read_index_speed                test_seqport
datatool                     readresult                      test_seqvector_ci
db95_bcp                     redoalignment_unit_test         test_serial
db95_cancel                  regexplocdemo                   test_server
db95_dbmorecmds              remote_app                      test_server_listeners
db95_done_handling           remote_app_client_sample        test-shift
db95_hang                    remote_blast_demo               test_snp_loader
db95_null                    remote_blast_unit_test          test_source_mod_parser
db95_null2                   remote_cgi                      test_sra
db95_numeric                 rmblastn                        test_sra_loader
db95_pending                 rpsblast                        test_stacktrace
db95_rpc                     rpstblastn                      test_staticmap
db95_setnull                 rps_unit_test                   test_strdbl
db95_t0001                   run_with_lock                   test_strsearch
db95_t0002                   scoreblk_unit_test              test_sub_reg
db95_t0003                   score_builder_unit_test         test_table_printer
db95_t0004                   sdbapi_advanced_features        test_tar
db95_t0005                   sdbapi_simple                   test-tds95
db95_t0006                   sdbapi_test_mt_pooling          test_tempstr
db95_t0007                   sdbapi_unit_test                test_threaded_client
db95_t0008                   search_strategy_unit_test       test_threaded_server
db95_t0009                   seedtop                         test_thread_pool
db95_t0011                   segmasker                       test_thread_pool_old
db95_t0012                   seqalign_unit_test              test_timsort
db95_t0013                   seqannot_splicer                test_title
db95_t0014                   seqdb_demo                      test_tls_object
db95_t0015                   seqdb_lmdb_unit_test            test_transmissionrw
db95_t0016                   seqdb_perf                      test_trial
db95_t0017                   seqdb_unit_test                 test_trial_check
db95_t0018                   seq_entry_reassign_ids          test_uncaught_exception
db95_t0019                   seqfeatdata_unit_test           test_uoconv
db95_t0020                   seqfeat_unit_test               test_user_agent
db95_t0021                   seq_id_unit_test                test_utf8
db95_t0022                   seqinfosrc_unit_test            test_uttp
db95_t0023                   seq_loc_unit_test               test_validator
db95_text_buffer             seqmasks_io_unit_test           test_value_convert
db95_thread                  seqsrc_unit_test                test_vdbgraph_loader
db95_timeout                 seqsub_split                    test_vmerge
dbapi_advanced_features      seqvec_bench                    test_weakref
dbapi_bcp                    setupfactory_unit_test          test_wgs_loader
dbapi_cache_admin            snp_test                        tls
dbapi_cache_test             soap_client_sample              tracebacksearch_unit_test
dbapi_conn_policy            soap_server_sample              traceback_unit_test
dbapi_context_test           socket_io_bouncer               uniform_search_unit_test
dbapi_cursor                 speedtest                       unit_test_5colftblreader
dbapi_driver_check           splign                          unit_test_agp_converter
dbapi_query                  split_cache                     unit_test_agp_seq_entry
dbapi_send_data              split_loader_demo               unit_test_align_cleanup
dbapi_simple                 split_query_unit_test           unit_test_alnmgr
dbapi_testspeed              sra_test                        unit_test_alnreader
dbapi_unit_test              srcchk                          unit_test_alnwriter
db_copy                      srsearch                        unit_test_alt_sample
deltablast                   stat_unit_test                  unit_test_autodef
delta_unit_test              streamtest                      unit_test_basic_cleanup
demo_contig_assembly         struct_dp_demo                  unit_test_bedgraphwriter
demo_gene_model              struct_util_demo                unit_test_bedreader
demo_genomic_compart         sub_cache_create                unit_test_bedwriter
demo_html                    subcheck                        unit_test_biosample_chk
demo_html_template           sub_fuse                        unit_test_biosample_util
demo_ncbi_clog               sub_image                       unit_test_bioseqgaps_ci
demo_score_builder           subj_ranges_unit_test           unit_test_cds_fix
demo_seqtest                 table2asn                       unit_test_defline
dump_asn_cache_index         tableval                        unit_test_entry_edit
dustmasker                   tblastn                         unit_test_eutils
ecnum_unit_test              tblastx                         unit_test_extended_cleanup
entrez2client                tds95_challenge                 unit_test_fasta_ostream
eutils_sample                tds95_charconv                  unit_test_fasta_reader
example_value_convert        tds95_collations                unit_test_feature_table_reader
fasthello.fcgi               tds95_condition                 unit_test_field
fcgi_sample.fcgi             tds95_convert                   unit_test_field_collection
feattree_sample              tds95_dataread                  unit_test_field_handler
feat_unit_test               tds95_dynamic1                  unit_test_format_guess
formatguess                  tds95_iconv_fread               unit_test_format_guess_ex
formatguess_unit_test        tds95_mutex1                    unit_test_gap_analysis
gapinfo_unit_test            tds95_nulls                     unit_test_gap_trim
gap_stats                    tds95_numeric                   unit_test_gene_model
gc_cli                       tds95_passarg                   unit_test_get_label
gencode_singleton_unit_test  tds95_readconf                  unit_test_gff3flybasewriter
gene_info_reader             tds95_strings                   unit_test_gff3reader
gene_info_unit_test          tds95_t0001                     unit_test_gff3reader_genbank
gene_info_writer             tds95_t0002                     unit_test_gff3writer
gene_info_writer_unit_test   tds95_t0003                     unit_test_gtfreader
genomic_compart_unit_test    tds95_t0004                     unit_test_gtfwriter
get_species_taxids.sh        tds95_t0005                     unit_test_gvfreader
gff_deconcat                 tds95_t0006                     unit_test_idmapper
gi2taxid                     tds95_t0007                     unit_test_id_mapper
graph_test                   tds95_t0008                     unit_test_internal_stops
grid_cgi_sample.cgi          tds95_toodynamic                unit_test_location_constraint
grid_cli                     tds95_utf8_1                    unit_test_loc_edit
grid_client_sample           tds95_utf8_2                    unit_test_mail_report
grid_worker_sample           tds95_utf8_3                    unit_test_mol_wt
gumbelparams                 test_algo_tree                  unit_test_objmgr
gumbelparams_unit_test       test_align                      unit_test_obj_sniff
hello.cgi                    test_annot_ci                   unit_test_orf
hfilter                      test_bam_loader                 unit_test_parse_text
hgvs2variation               test_base64                     unit_test_pmcidconv_client
hooks_commented              test_bdb_cursor                 unit_test_polya
hooks_copy_member            test_biotree                    unit_test_pub_edit
hooks_copy_object            test_bm                         unit_test_rna_edit
hooks_copy_variant           test_buffer_writer              unit_test_sample
hooks_highest_se_objs        test_bulkinfo                   unit_test_score_builder
hooks_read_member            test_bulkinfo_mt                unit_test_seq
hooks_read_object            test_cache_mt                   unit_test_seq_loc_util
hooks_read_variant           test_cgi_entry_reader           unit_test_seq_translator
hooks_skip_member            test_chainer                    unit_test_seq_trimmer
hooks_skip_object            test_checksum                   unit_test_so_map
hooks_skip_variant           test_compound_id                unit_test_source_edit
hooks_write_member           test_compress                   unit_test_splign
hooks_write_object           test_compress_archive           unit_test_srcwriter
hooks_write_variant          test_compress_mt                unit_test_string_constraint
hspfilter_besthit_unit_test  test_condvar                    unit_test_ucscreader
hspfilter_culling_unit_test  test_conn_stream_pushback       unit_test_ucscwriter
hspstream_unit_test          test_conn_tar                   unit_test_validator
http_connector_hit           test_csra_loader                unit_test_vcfreader
http_session_sample          test_csra_loader_mt             unit_test_vcfwriter
id1_fetch                    test-ct95                       unit_test_wigreader
id1_fetch_simple             test_ctransition                unit_test_wigwriter
id2_fetch_simple             test_ctransition_nlmzip         update_blastdb.pl
idfetch                      test_date                       update_prot_id
idmapper                     test-db95                       vdb_test
id_unit_test                 test_diag_parser                vecscreen
id_unit_test_bad             test_diff                       version_reference_unit_test
igblastn                     test_edit_saver                 vsrun_sample
igblastp                     test_eutils_client              walk_cache_test
image_info                   test_expr                       wgs_test
lang_query                   test_fasta_round_trip           wig2table
lds2_indexer                 test_feat_overlap               windowmasker
lds2_sample                  test_feat_tree                  windowmasker_2.2.22_adapter.py
legacy_blast.pl              test_floating_point_comparison  writedb_lmdb_unit_test
linkhsp_unit_test            test_fstream_pushback           writedb_unit_test
lmdb_test1                   test_fw                         xcompareannotsdemo
lmdb_test2                   test_get_console_password
lmdb_test3                   test_gridclient_stress

```







