

document.querySelector('title').textContent = 'SMRT Analysis on Biowulf';
SMRT Analysis on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Interactive job on Biowulf](#int)
[Batch job on Biowulf](#sbatch)
 |



SMRT® Analysis is a bioinformatics software suite available for analysis of DNA sequencing data from Pacific Biosciences’ SMRT technology. Users can choose from a variety of analysis protocols that utilize PacBio® and third-party tools. Analysis protocols include de novo genome assembly, cDNA mapping, DNA base-modification detection, and long-amplicon analysis to determine phased consensus sequences.



Documentation
* [SMRT Link main site](https://www.pacb.com/support/software-downloads/) (SMRT Analysis is part of SMRT Link)
+ See, in particular, the SMRT tools reference guide and the SMRT Link user guide linked from there.

* [SMRT Link on GitHub](https://github.com/PacificBiosciences/SMRT-Link)


Important Notes
* Module Name: smrtanalysis (see [the modules page](/apps/modules.html) for more information)
* Environment variables set 
	+ SMRT\_HOME* Barcodes in /fdb/smrtanalysis/barcodes* Sample data in /fdb/smrtanalysis/canneddata



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.

This is a sample interactive session of the lambda phage site acceptance test done on the local node using pbcromwell. (user input in **bold**):



```

[teacher@biowulf ~]$ **sinteractive --cpus-per-task=12**
salloc.exe: Pending job allocation 43027948
salloc.exe: job 43027948 queued and waiting for resources
salloc.exe: job 43027948 has been allocated resources
salloc.exe: Granted job allocation 43027948
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3109 are ready for job
srun: error: x11: no local DISPLAY defined, skipping
[teacher@cn3109 smrtanalysis]$ **module load smrtanalysis**
[+] Loading smrtanalysis 8.0.0.79519
[teacher@cn3109 ~]$ **mkdir /data/$USER/smrtanalysis**
[teacher@cn3109 ~]$ **cd !$**
[teacher@cn3109 smrtanalysis]$ **pbcromwell configure**
[INFO] 2019-10-30 18:32:18,797Z [pbcromwell.cli] Using pbcommand v1.9.2
[INFO] 2019-10-30 18:32:18,798Z [pbcromwell.cli] completed setting up logger with 
[INFO] 2019-10-30 18:32:18,798Z [pbcromwell.cli] log opts {'file\_name': None, 'level': 20}
[WARNING] 2019-10-30 18:32:18,798Z [pbcromwell.cli] No database port specified - will run with in-memory DB
[INFO] 2019-10-30 18:32:18,800Z [pbcromwell.cli] Wrote config file to cromwell.conf
[INFO] 2019-10-30 18:32:18,801Z [pbcromwell.cli] exiting with return code 0 in 0.00 sec.
[teacher@cn3109 smrtanalysis]$ **ls**
cromwell.conf
[teacher@cn3109 smrtanalysis]$ **pbcromwell show-workflows**
[INFO] 2019-10-30 18:33:09,770Z [pbcromwell.cli] Using pbcommand v1.9.2
[INFO] 2019-10-30 18:33:09,770Z [pbcromwell.cli] completed setting up logger with 
[INFO] 2019-10-30 18:33:09,770Z [pbcromwell.cli] log opts {'file\_name': None, 'level': 20}
[INFO] 2019-10-30 18:33:09,770Z [pbcromwell.cli] SMRT\_PIPELINE\_BUNDLE\_DIR="/usr/local/apps/smrtanalysis/8.0.0.79519/install/smrtlink-release\_8.0.0.79519/bundles/smrttools/install/smrttools-release\_8.0.0.79510/private/thirdparty/python/python\_2.7.16/binwrap/../../../../pacbio/pbpipeline-resources"


cromwell.workflows.pb\_hgap4: Assembly (HGAP4)
cromwell.workflows.pb\_basemods: Base Modification Analysis
cromwell.workflows.pb\_ccs\_mapping: CCS with Mapping
cromwell.workflows.pb\_ccs: Circular Consensus Sequencing (CCS)
cromwell.workflows.pb\_bam2fastx: Convert BAM to FASTX
cromwell.workflows.pb\_demux\_ccs: Demultiplex Barcodes
cromwell.workflows.pb\_demux\_subreads: Demultiplex Barcodes
cromwell.workflows.pb\_isoseq3: Iso-Seq
cromwell.workflows.pb\_isoseq3\_ccsonly: Iso-Seq
cromwell.workflows.pb\_laa: Long Amplicon Analysis (LAA)
cromwell.workflows.pb\_align\_ccs: Mapping
cromwell.workflows.pb\_assembly\_microbial: Microbial Assembly
cromwell.workflows.pb\_mv\_ccs: Minor Variants Analysis
cromwell.workflows.pb\_resequencing: Resequencing
cromwell.workflows.pb\_sat: Site Acceptance Test (SAT)
cromwell.workflows.pb\_sv\_ccs: Structural Variant Calling
cromwell.workflows.pb\_sv\_clr: Structural Variant Calling

Run 'pbcromwell show-workflow-details ' to display further
information about a workflow. Note that the cromwell.workflows.
prefix is optional.
[INFO] 2019-10-30 18:33:09,840Z [pbcromwell.cli] exiting with return code 0 in 0.07 sec.
[teacher@cn3109 smrtanalysis]$ **pbcromwell show-workflow-details pb\_sat**
[INFO] 2019-10-30 18:33:36,290Z [pbcromwell.cli] Using pbcommand v1.9.2
[INFO] 2019-10-30 18:33:36,290Z [pbcromwell.cli] completed setting up logger with 
[INFO] 2019-10-30 18:33:36,290Z [pbcromwell.cli] log opts {'file\_name': None, 'level': 20}
[INFO] 2019-10-30 18:33:36,290Z [pbcromwell.cli] SMRT\_PIPELINE\_BUNDLE\_DIR="/usr/local/apps/smrtanalysis/8.0.0.79519/install/smrtlink-release\_8.0.0.79519/bundles/smrttools/install/smrttools-release\_8.0.0.79510/private/thirdparty/python/python\_2.7.16/binwrap/../../../../pacbio/pbpipeline-resources"


Pipeline Summary
Pipeline Id: cromwell.workflows.pb\_sat
Name : Site Acceptance Test (SAT)
Description: Cromwell workflow pb\_sat
EntryPoints: 2
 eid\_ref\_dataset -> PacBio.DataSet.ReferenceSet
 eid\_subread -> PacBio.DataSet.SubreadSet
Tags : cromwell, mapping
Task Options:
 dataset\_filters =
 downsample\_factor = 0

[INFO] 2019-10-30 18:33:36,291Z [pbcromwell.cli] exiting with return code 0 in 0.00 sec.
[teacher@cn3109 smrtanalysis]$ **pbcromwell run pb\_sat \
 --entry /fdb/smrtanalysis/canneddata/lambdaTINY/m54026\_181219\_010936\_tiny.subreadset.xml \
 --entry /fdb/smrtanalysis/canneddata/referenceset/lambdaNEB/referenceset.xml \
 --config $PWD/cromwell.conf \
 --nproc $SLURM\_CPUS\_PER\_TASK \
 --output-dir sat** 

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).

While pbcromwell can be run with the slurm backend to submit jobs while the supervisor process runs in an interactive session, we recommend queuing the supervisor process itself. Create a batch input file (e.g. smrtanalysis.sh) as follows:




```

#!/bin/bash
set -e

module load smrtanalysis

### cromwell configuration
pbcromwell configure
# The generated cromwell.conf cannot be used as is. The following command makes the necessary adjustments.
sed -i -r \
  -e 's%(sbatch|scancel|squeue)%/usr/local/slurm/bin/\1%' `# use absolute paths for slurm commands (slurm directory is purged from $PATH by PacBio wrapper)` \
  -e '/job-id-regex/s/(job-id-regex\s+=\s+).*/\1 "(\\\\d+)"/' `# NIH HPC's sbatch has non-standard output` \
  cromwell.conf

### run workflow
pbcromwell run pb_sat \
 --entry /fdb/smrtanalysis/canneddata/lambdaTINY/m54026_181219_010936_tiny.subreadset.xml \
 --entry /fdb/smrtanalysis/canneddata/referenceset/lambdaNEB/referenceset.xml \
 --config $PWD/cromwell.conf \
 --backend slurm \
 --output-dir sat 

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command, requesting a walltime of the expected run length of the whole pipeline.



```
sbatch [--time=#] smrtanalysis.sh
```







