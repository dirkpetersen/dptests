

document.querySelector('title').textContent = 'FCS on Biowulf';
FCS on Biowulf


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


 FCS is a toolset to elimate contaminant sequences from a genome assembly. Currently, there are two tools in this toolset: `fcs-adaptor` and `fcs-gx`. 


### References:


* Astashyn A, Tvedte ES, Sweeney D, Sapojnikov V, Bouk N, Joukov V, Mozes E, Strope PK, Sylla PM, Wagner L, Bidwell SL, Clark K, Davis EW, Smith-White B, Hlavina W, Pruitt KD, Schneider VA, Murphy TD
*[Rapid and sensitive detection of genome contamination at scale with FCS-GX](https://pubmed.ncbi.nlm.nih.gov/37292984/)*. bioRxiv [Preprint]. 2023 Jun 6:2023.06.02.543519. doi: 10.1101/2023.06.02.543519.


Documentation
* FCS Main Site: [Github site](https://github.com/ncbi/fcs)


Important Notes
* Module Name: fcs (see [the modules page](/apps/modules.html) for more information)
* Test data can be found in `${FCS_TEST_DATA}`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. From the FCS tutorials at the [FCS Github site](https://github.com/ncbi/fcs):



```

[user@biowulf]$ **sinteractive --mem=3g --cpus-per-task=4 --gres=lscratch:100**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn4224 are ready for job

[user@cn4224 ~]$ **module load fcs**
[+] Loading fcs  0.4.0  on cn4285
[+] Loading singularity  3.10.5  on cn4285

[user@cn4244 ~]$ **cd /lscratch/${SLURM\_JOB\_ID}**
[user@cn4224 ~]$ **cp ${FCS\_TEST\_DATA}/fcsadaptor\_prok\_test.fa.gz inputdir/.**
[user@cn4224 ~]$ **mkdir inputdir outputdir**
[user@cn4224 ~]$ **run\_fcsadaptor.sh --fasta-input ./inputdir/fcsadaptor\_prok\_test.fa.gz --output-dir ./outputdir --prok --container-engine singularity --image ${FCS\_HOME}/fcs-adaptor.sif**
[WARN  tini (2647065)] Tini is not running as PID 1 and isn't registered as a child subreaper.
Zombie processes will not be re-parented to Tini, so zombie reaping won't work.
To fix the problem, use the -s option or set the environment variable TINI_SUBREAPER to register Tini as a child subreaper, or run Tini as PID 1.
Output will be placed in: /output-volume
Executing the workflow
Resolved '/app/fcs/progs/ForeignContaminationScreening.cwl' to 'file:///app/fcs/progs/ForeignContaminationScreening.cwl'
[workflow ] start
[workflow ] starting step ValidateInputSequences
[step ValidateInputSequences] start

[..]

[job all_skipped_trims] completed success
[step all_skipped_trims] completed success
[workflow ] starting step all_cleaned_fasta
[step all_cleaned_fasta] start
[step all_cleaned_fasta] completed success
[workflow ] completed success

[user@cn4224 ~]$ **cp ${FCS\_TEST\_DATA}/fcsgx\_test.fa.gz .**
[user@cn4224 ~]$ **SOURCE\_DB\_MANIFEST="https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.manifest"**
[user@cn4224 ~]$ **LOCAL\_DB=/lscratch/${SLURM\_JOB\_ID}/gxdb**
[user@cn4224 ~]$ **fcs.py db get --mft "$SOURCE\_DB\_MANIFEST" --dir "$LOCAL\_DB/test-only"**
[user@cn4224 ~]$ **fcs.py db check --mft "$SOURCE\_DB\_MANIFEST" --dir "$LOCAL\_DB/gxdb"**
===============================================================================
Source:      https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only
Destination: /app/db/gxdb
Resuming failed transfer in /app/db/gxdb...
Space check: Available:15.71GiB; Existing:0B; Incoming:4.29GiB; Delta:4.29GiB

Requires transfer: 56B test-only.meta.jsonl

Requires transfer: 5.92KiB test-only.taxa.tsv

Requires transfer: 21.33KiB test-only.seq_info.tsv.gz

Requires transfer: 7.85MiB test-only.blast_div.tsv.gz

Requires transfer: 67.56MiB test-only.gxs

Requires transfer: 4.21GiB test-only.gxi

[user@cn4224 ~]$ **GXDB\_LOC=/lscratch/${SLURM\_JOB\_ID}/gxdb**
[user@cn4224 ~]$ **fcs.py screen genome --fasta ./fcsgx\_test.fa.gz --out-dir ./gx\_out/ --gx-db "$GXDB\_LOC/test-only" --tax-id 6973**
--------------------------------------------------------------------

tax-id    : 6973
fasta     : /sample-volume/fcsgx_test.fa.gz
size      : 8.55 MiB
split-fa  : True
BLAST-div : roaches
gx-div    : anml:insects
w/same-tax: True
bin-dir   : /app/bin
gx-db     : /app/db/gxdb/test-only/test-only.gxi
gx-ver    : Mar 10 2023 15:34:33; git:v0.4.0-3-g8096f62
output    : /output-volume//fcsgx_test.fa.6973.taxonomy.rpt

--------------------------------------------------------------------

 [...]

fcs_gx_report.txt contamination summary:
----------------------------------------
                                seqs      bases
                               ----- ----------
TOTAL                            243   27170378
-----                          ----- ----------
prok:CFB group bacteria          243   27170378

--------------------------------------------------------------------

fcs_gx_report.txt action summary:
---------------------------------
                                seqs      bases
                               ----- ----------
TOTAL                            243   27170378
-----                          ----- ----------
EXCLUDE                          214   25795430
REVIEW                            29    1374948

--------------------------------------------------------------------

[user@cn4224 ~]$ **zcat fcsgx\_test.fa.gz | fcs.py clean genome --action-report ./gx\_out/fcsgx\_test.fa.6973.fcs\_gx\_report.txt --output clean.fasta --contam-fasta-out contam.fasta**
Applied 214 actions; 25795430 bps dropped; 0 bps hardmasked.

[user@cn4224 ~]$ **ls gx\_out**
fcsgx_test.fa.6973.fcs_gx_report.txt  fcsgx_test.fa.6973.taxonomy.rpt


```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. fcsadaptor.sh) similar to the following.



```

#! /bin/bash

module load fcs
cd /lscratch/${SLURM_JOB_ID}
mkdir inputdir outputdir
cp ${FCS_TEST_DATA}/fcsadaptor_prok_test.fa.gz inputdir/.

run_fcsadaptor.sh --fasta-input ./inputdir/fcsadaptor_prok_test.fa.gz --output-dir ./outputdir --prok --container-engine singularity --image ${FCS_HOME}/fcs-adaptor.sif

cp -r outputdir /data/${USER}/.

```

Submit these jobs using the Slurm [sbatch](/docs/userguide.html) command.


Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile to run fcs-adaptor (e.g. fcsadaptor.swarm). For example:



```

run_fcsadaptor.sh --fasta-input ./inputdir/fcsadaptor_prok_test1.fa.gz --output-dir ./outputdir1 --prok --container-engine singularity --image ${FCS_HOME}/fcs-adaptor.sif; cp -r outputdir1 /data/${USER}/. 
run_fcsadaptor.sh --fasta-input ./inputdir/fcsadaptor_prok_test2.fa.gz --output-dir ./outputdir2 --prok --container-engine singularity --image ${FCS_HOME}/fcs-adaptor.sif; cp -r outputdir2 /data/${USER}/.
run_fcsadaptor.sh --fasta-input ./inputdir/fcsadaptor_prok_test3.fa.gz --output-dir ./outputdir3 --prok --container-engine singularity --image ${FCS_HOME}/fcs-adaptor.sif; cp -r outputdir3 /data/${USER}/.

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f fcsadaptor.swarm [-g #] --module fcs
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --module fcs  Loads the fcs module for each subjob in the swarm
 | |
 | |








