

document.querySelector('title').textContent = 'AnnotSV on Biowulf';
AnnotSV on Biowulf


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



AnnotSV is a program designed for annotating Structural Variations (SV). This tool compiles functionally, regulatory and clinically relevant information and aims at providing annotations useful to i) interpret SV potential pathogenicity and ii) filter out SV potential false positives.



Documentation
* [AnnotSV Main Site](http://lbgi.fr/AnnotSV/)


Important Notes
* Module Name: annotsv (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded app
* environment variables set 
	+ ANNOTSV* Example files in /usr/local/apps/annotsv/Example



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf ~]$ **sinteractive -c4 --mem=10g --gres=lscratch:10**
salloc.exe: Pending job allocation 64198846
salloc.exe: job 64198846 queued and waiting for resources
salloc.exe: job 64198846 has been allocated resources
salloc.exe: Granted job allocation 64198846
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0890 are ready for job
srun: error: x11: no local DISPLAY defined, skipping

[user@cn0890 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn0890 64198846]$ **module load annotsv**
[+] Loading annotsv  2.2  on cn0890
[+] Loading tcl_tk 8.6.8
[+] Loading bedtools  2.29.2

[user@cn0890 64198846]$ **cp ${ANNOTSV}/Example/\* .**

[user@cn0890 64198846]$ **ls**
commands.README  test.bed  test_GRCh37.annotated.example.tsv

[user@cn0890 64198846]$ **AnnotSV -SVinputFile test.bed -SVinputInfo 1 -outputFile ./test.annotated.tsv**
AnnotSV 2.2

Copyright (C) 2017-2019 GEOFFROY Veronique

Please feel free to contact me for any suggestions or bug reports
email: veronique.geoffroy@inserm.fr

Tcl/Tk version: 8.6

[snip...]

...annotation in progress (September 03 2020 - 12:13)


...Output columns annotation:
        AnnotSV ID; SV chrom; SV start; SV end; SV length; ; ; AnnotSV type; Gene name; NM; CDS length; tx length; location; location2; intersectStart; intersectEnd; DGV_GAIN_IDs; DGV_GAIN_n_samples_with_SV; DGV_GAIN_n_samples_tested; DGV_GAIN_Frequency; DGV_LOSS_IDs; DGV_LOSS_n_samples_with_SV; DGV_LOSS_n_samples_tested; DGV_LOSS_Frequency; GD_ID; GD_AN; GD_N_HET; GD_N_HOMALT; GD_AF; GD_POPMAX_AF; GD_ID_others; DDD_SV; DDD_DUP_n_samples_with_SV; DDD_DUP_Frequency; DDD_DEL_n_samples_with_SV; DDD_DEL_Frequency; 1000g_event; 1000g_AF; 1000g_max_AF; IMH_ID; IMH_AF; IMH_ID_others; promoters; dbVar_event; dbVar_variant; dbVar_status; TADcoordinates; ENCODEexperiments; GCcontent_left; GCcontent_right; Repeats_coord_left; Repeats_type_left; Repeats_coord_right; Repeats_type_right; ACMG; HI_CGscore; TriS_CGscore; DDD_status; DDD_mode; DDD_consequence; DDD_disease; DDD_pmids; HI_DDDpercent; synZ_ExAC; misZ_ExAC; pLI_ExAC; delZ_ExAC; dupZ_ExAC; cnvZ_ExAC; morbidGenes; morbidGenesCandidates; Mim Number; Phenotypes; Inheritance


...AnnotSV is done with the analysis (September 03 2020 - 12:13)

[user@cn0890 64198846]$ **exit**
exit
salloc.exe: Relinquishing job allocation 64198846
salloc.exe: Job allocation 64198846 has been revoked.

[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. annotsv.sh). For example:



```

#!/bin/bash
module load annotsv
AnnotSV -SVinputFile test.bed -SVinputInfo 1 -outputFile ./test.annotated.tsv

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch annotsv.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. annotsv.swarm). For example:



```

AnnotSV -SVinputFile test1.bed -SVinputInfo 1 -outputFile ./test1.annotated.tsv
AnnotSV -SVinputFile test2.bed -SVinputInfo 1 -outputFile ./test2.annotated.tsv
AnnotSV -SVinputFile test3.bed -SVinputInfo 1 -outputFile ./test3.annotated.tsv

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f annotsv.swarm -g 10 --module annotsv
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module annotsv Loads the annotsv module for each subjob in the swarm 
 | |
 | |
 | |








