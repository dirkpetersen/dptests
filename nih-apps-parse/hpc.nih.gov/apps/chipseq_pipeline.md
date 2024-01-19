

document.querySelector('title').textContent = 'chipseq\_pipeline on Biowulf';
chipseq\_pipeline on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
 |



This ChIP-Seq pipeline is based off the ENCODE (phase-3) transcription factor and histone ChIP-seq pipeline specifications.



Documentation
* [Main GitHub repository](https://github.com/ENCODE-DCC/chip-seq-pipeline2)


Important Notes
* Module Name: chipseq\_pipeline (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded app (using WDL)
 * Reference data in /fdb/chipseq\_pipeline/<version>* Environment variables set 
	
	 Version 1.2.0
	 
		- CROMWELL\_JAR=/usr/local/apps/cromwell/40/cromwell-40.jar
		- CHIPSEQ\_HOME
		- CHIPSEQ\_TESTDATA
	 Versions 1.4.0.1 / 1.5.1 / 1.6.1 / 2.1.6
	 
		- CSP\_WDL
		- CSP\_TEST\_JAR* Version 2.1.6+ uses caper 2.1 which has requires new config. You should backup and replace your ~/.caper/ configuration prior to use.
 * Version 2.1.6+ now requires the --singularity flag in caper run matching other ENCODE pipelines available on Biowulf.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



* 1.2.0
* 1.4.0.1 / 1.5.1 / 1.6.1
* 2.1.6





```


[user@biowulf]$ **sinteractive -c 12 --mem 20g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load chipseq\_pipeline**

[user@cn3144 ~]$ **cp -r $CHIPSEQ\_TESTDATA/\* .**

[user@cn3144 ~]$ **INPUT=ENCSR936XTK\_subsampled\_chr19\_only.json**

[user@cn3144 ~]$ **java -jar -Xmx15G -Dconfig.file=$CHIPSEQ\_HOME/backends/backend.conf
-Dbackend.default=singularity $CROMWELL\_JAR run $CHIPSEQ\_HOME/chip.wdl
-i ${INPUT} -o $CHIPSEQ\_HOME/workflow\_opts/singularity.json -m meta.json**
[2019-06-04 12:22:50,83] [info] Running with database db.url = jdbc:hsqldb:mem:d0de5da6-b561-40fe-8cb9-a1b3cf4edb0b;shutdown=false;hsqldb.tx=mvcc
[2019-06-04 12:22:58,86] [info] Running migration RenameWorkflowOptionsInMetadata with a read batch size of 100000 and a write batch size of 100000
[2019-06-04 12:22:58,88] [info] [RenameWorkflowOptionsInMetadata] 100%
[2019-06-04 12:22:59,01] [info] Running with database db.url = jdbc:hsqldb:mem:e0b3ccbd-65dd-41da-923d-2f806691bbab;shutdown=false;hsqldb.tx=mvcc
[2019-06-04 12:22:59,41] [warn] This actor factory is deprecated. Please use cromwell.backend.google.pipelines.v1alpha2.PipelinesApiLifecycleActorFactory for PAPI v1 or cromwell.backend.google.pipelines.v2alpha1.PipelinesApiLifecycleActorFactory for PAPI v2
[2019-06-04 12:22:59,69] [info] Slf4jLogger started
[2019-06-04 12:22:59,97] [info] Workflow heartbeat configuration:
{
  "cromwellId" : "cromid-3b31626",
  "heartbeatInterval" : "2 minutes",
  "ttl" : "10 minutes",
  "writeBatchSize" : 10000,
  "writeThreshold" : 10000
}
[2019-06-04 12:23:00,07] [info] Metadata summary refreshing every 1 second.
[2019-06-04 12:23:00,11] [info] KvWriteActor configured to flush with batch size 200 and process rate 5 seconds.
[2019-06-04 12:23:00,12] [info] WriteMetadataActor configured to flush with batch size 200 and process rate 5 seconds.
[2019-06-04 12:23:00,12] [info] CallCacheWriteActor configured to flush with batch size 100 and process rate 3 seconds.
[2019-06-04 12:23:00,12] [warn] 'docker.hash-lookup.gcr-api-queries-per-100-seconds' is being deprecated, use 'docker.hash-lookup.gcr.throttle' instead (see reference.conf)
[...]
[2019-06-04 13:12:45,82] [info] Connection pools shut down
[2019-06-04 13:12:45,82] [info] Shutting down SubWorkflowStoreActor - Timeout = 1800 seconds
[2019-06-04 13:12:45,82] [info] Shutting down JobStoreActor - Timeout = 1800 seconds
[2019-06-04 13:12:45,82] [info] Shutting down CallCacheWriteActor - Timeout = 1800 seconds
[2019-06-04 13:12:45,82] [info] Shutting down ServiceRegistryActor - Timeout = 1800 seconds
[2019-06-04 13:12:45,82] [info] Shutting down DockerHashActor - Timeout = 1800 seconds
[2019-06-04 13:12:45,82] [info] Shutting down IoProxy - Timeout = 1800 seconds
[2019-06-04 13:12:45,82] [info] SubWorkflowStoreActor stopped
[2019-06-04 13:12:45,82] [info] CallCacheWriteActor Shutting down: 0 queued messages to process
[2019-06-04 13:12:45,82] [info] CallCacheWriteActor stopped
[2019-06-04 13:12:45,82] [info] WriteMetadataActor Shutting down: 0 queued messages to process
[2019-06-04 13:12:45,82] [info] KvWriteActor Shutting down: 0 queued messages to process
[2019-06-04 13:12:45,82] [info] JobStoreActor stopped
[2019-06-04 13:12:45,83] [info] IoProxy stopped
[2019-06-04 13:12:45,83] [info] ServiceRegistryActor stopped
[2019-06-04 13:12:45,83] [info] DockerHashActor stopped
[2019-06-04 13:12:45,84] [info] Database closed
[2019-06-04 13:12:45,84] [info] Stream materializer shut down
[2019-06-04 13:12:45,84] [info] WDL HTTP import resolver closed


[user@cn3144 ~]$ **tail meta.json**
    "chip.idr_thresh": 0.05
  },
  "labels": {
    "cromwell-workflow-id": "cromwell-b3a17a9a-08aa-491b-ba91-cffaf0346068"
  },
  "submission": "2019-06-04T13:34:43.516-04:00",
  "status": "Succeeded",
  "end": "2019-06-04T14:11:36.841-04:00",
  "start": "2019-06-04T13:34:43.589-04:00"
}[apptest4@cn3585 test]$ 
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





```


[user@biowulf ~]$ **sinteractive --cpus-per-task=12 --mem=15g --gres=lscratch:10**
salloc.exe: Pending job allocation 56941209
salloc.exe: job 56941209 queued and waiting for resources
salloc.exe: job 56941209 has been allocated resources
salloc.exe: Granted job allocation 56941209
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3130 are ready for job
srun: error: x11: no local DISPLAY defined, skipping

[user@cn3130 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn3130 56941209]$ **module load chipseq\_pipeline/<ver>**
[+] Loading chipseq_pipeline  <ver>  on cn3130

[user@cn3130 56941209]$ **cp $CSP\_TEST\_JAR .**

[user@cn3130 56941209]$ **mkdir -p ~/.caper**

[user@cn3130 56941209]$ **caper init local**
Downloading Cromwell JAR... https://github.com/broadinstitute/cromwell/releases/download/47/cromwell-47.jar
Downloading Womtool JAR... https://github.com/broadinstitute/cromwell/releases/download/47/womtool-47.jar

[user@cn3130 56941209]$ **caper run $CSP\_WDL -i ENCSR000DYI\_subsampled\_chr19\_only.json**
2020-04-27 13:47:09,505|autouri.autouri|INFO| cp: copied, src=https://storage.googleapis.com/encode-pipeline-genome-data/hg38_chr19_chrM/GRCh38_no_alt_analysis_set_GCA_000001405.15.chr19_chrM.fasta.gz, dest=/lscratch/56941209/.caper_tmp/dcc664bb3e886aa9127f3bd8522bff3c/GRCh38_no_alt_analysis_set_GCA_000001405.15.chr19_chrM.fasta.gz
[...snip...]
2020-04-27 14:07:55,306 cromwell-system-akka.dispatchers.engine-dispatcher-43 INFO  - SingleWorkflowRunnerActor workflow finished with status 'Succeeded'.
{
  "outputs": {
    "chip.report": "/lscratch/56941209/chip/e4dac084-f786-4f24-a567-4c63e34a1040/call-qc_report/execution/glob-eae855c82d0f7e2185388856e7b2cc7b/qc.html",
    "chip.qc_json_ref_match": false,
    "chip.qc_json": "/lscratch/56941209/chip/e4dac084-f786-4f24-a567-4c63e34a1040/call-qc_report/execution/glob-3440f922973abb7a616aaf203e0db08b/qc.json"
  },
  "id": "e4dac084-f786-4f24-a567-4c63e34a1040"
}
2020-04-27 14:07:58,575 cromwell-system-akka.dispatchers.engine-dispatcher-103 INFO  - SingleWorkflowRunnerActor writing metadata to /lscratch/56941209/.caper_tmp/chip/20200427_134708_801102/metadata.json
2020-04-27 14:07:58,631  INFO  - Workflow polling stopped
2020-04-27 14:07:58,648  INFO  - 0 workflows released by cromid-c5dda4a
2020-04-27 14:07:58,651  INFO  - Shutting down WorkflowStoreActor - Timeout = 5 seconds
2020-04-27 14:07:58,655  INFO  - Shutting down WorkflowLogCopyRouter - Timeout = 5 seconds
2020-04-27 14:07:58,658  INFO  - Shutting down JobExecutionTokenDispenser - Timeout = 5 seconds
2020-04-27 14:07:58,662  INFO  - JobExecutionTokenDispenser stopped
2020-04-27 14:07:58,662 cromwell-system-akka.dispatchers.engine-dispatcher-104 INFO  - Aborting all running workflows.
2020-04-27 14:07:58,663  INFO  - WorkflowStoreActor stopped
2020-04-27 14:07:58,671  INFO  - WorkflowLogCopyRouter stopped
2020-04-27 14:07:58,671  INFO  - Shutting down WorkflowManagerActor - Timeout = 3600 seconds
2020-04-27 14:07:58,672 cromwell-system-akka.dispatchers.engine-dispatcher-80 INFO  - WorkflowManagerActor All workflows finished
2020-04-27 14:07:58,672  INFO  - WorkflowManagerActor stopped
2020-04-27 14:07:59,126  INFO  - Connection pools shut down
2020-04-27 14:07:59,132  INFO  - Shutting down SubWorkflowStoreActor - Timeout = 1800 seconds
2020-04-27 14:07:59,132  INFO  - Shutting down JobStoreActor - Timeout = 1800 seconds
2020-04-27 14:07:59,132  INFO  - Shutting down CallCacheWriteActor - Timeout = 1800 seconds
2020-04-27 14:07:59,132  INFO  - Shutting down ServiceRegistryActor - Timeout = 1800 seconds
2020-04-27 14:07:59,132  INFO  - Shutting down DockerHashActor - Timeout = 1800 seconds
2020-04-27 14:07:59,132  INFO  - Shutting down IoProxy - Timeout = 1800 seconds
2020-04-27 14:07:59,132  INFO  - SubWorkflowStoreActor stopped
2020-04-27 14:07:59,134  INFO  - CallCacheWriteActor Shutting down: 0 queued messages to process
2020-04-27 14:07:59,134  INFO  - WriteMetadataActor Shutting down: 0 queued messages to process
2020-04-27 14:07:59,134  INFO  - CallCacheWriteActor stopped
2020-04-27 14:07:59,134  INFO  - KvWriteActor Shutting down: 0 queued messages to process
2020-04-27 14:07:59,142  INFO  - JobStoreActor stopped
2020-04-27 14:07:59,142  INFO  - IoProxy stopped
2020-04-27 14:07:59,142  INFO  - ServiceRegistryActor stopped
2020-04-27 14:07:59,143  INFO  - DockerHashActor stopped
2020-04-27 14:07:59,167  INFO  - Database closed
2020-04-27 14:07:59,168  INFO  - Stream materializer shut down
2020-04-27 14:07:59,170  INFO  - WDL HTTP import resolver closed
2020-04-27 14:07:59,434|caper.caper|INFO| run: 0, e4dac084-f786-4f24-a567-4c63e34a1040, None

[user@cn3130 56941209]$ **exit**
exit
salloc.exe: Relinquishing job allocation 56941209
salloc.exe: Job allocation 56941209 has been revoked.

```





```


[user@biowulf ~]$ **sinteractive --cpus-per-task=12 --mem=15g --gres=lscratch:30**
salloc: Pending job allocation 41186334
salloc: job 41186334 queued and waiting for resources
salloc: job 41186334 has been allocated resources
salloc: Granted job allocation 41186334
salloc: Waiting for resource configuration
salloc: Nodes cn0873 are ready for job
srun: error: x11: no local DISPLAY defined, skipping
error: unable to open file /tmp/slurm-spank-x11.41186334.0
slurmstepd: error: x11: unable to read DISPLAY value

[user@cn0873 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn0873 41186334]$ **module load chipseq\_pipeline/<ver>**
[+] Loading chipseq_pipeline  <ver>  on cn0873

        
```

If you have used caper for an earlier version of an ENCODE pipeline you should regenerate your caper config. In this example the pipeline will only be run locally - i.e. it will not submit
 tasks as slurm jobs. Follow the [caper](https://github.com/ENCODE-DCC/caper)
 docs to set up a config file for slurm submission. This has to be done only once.



```

[user@cn0873 41186334]$ **[[ -d ~/.caper ]] && mv ~/.caper ~/caper.$(date +%F).bak # back up old caper config**        
        
[user@cn0873 41186334]$ **mkdir -p ~/.caper && caper init local**
        
```


```

[user@cn0873 41186334]$ **caper run $CSP\_WDL -i $CSP\_TEST\_JAR --singularity**
2022-06-06 14:15:50,807|caper.cli|INFO| Cromwell stdout: /lscratch/41186334/cromwell.out
2022-06-06 14:15:50,807|caper.caper_base|INFO| Creating a timestamped temporary directory. /lscratch/41186334/.caper_tmp/chip/20220606_135150_807598
2022-06-06 14:15:50,807|caper.caper_runner|INFO| Localizing files on work_dir. /lscratch/41186334/.caper_tmp/chip/20220606_135150_807598
[...snip...]
2022-06-06 14:43:59,733|caper.cromwell_workflow_monitor|INFO| Workflow: id=dcee6b44-b2bd-414e-897b-baddc2d746b2, status=Succeeded
2022-06-06 14:44:26,941|caper.cromwell_metadata|INFO| Wrote metadata file. /lscratch/41186334/chip/dcee6b44-b2bd-414e-897b-baddc2d746b2/metadata.json
2022-06-06 14:44:26,942|caper.nb_subproc_thread|INFO| Cromwell finished successfully.

[user@cn0873 41186334]$ **exit**
exit
salloc.exe: Relinquishing job allocation 41186334
salloc.exe: Job allocation 41186334 has been revoked.

```




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. chipseq\_pipeline.sh). For example:




```

#!/bin/bash
set -e
module load chipseq_pipeline
cp $CSP_TEST_JAR .
INPUT=ENCSR000DYI_subsampled_chr19_only.json
caper run $CSP_WDL -i $INPUT # version 2+ requires --singularity
        
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=12 --mem=15g chipseq_pipeline.sh
```







