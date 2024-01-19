

document.querySelector('title').textContent = 'MEGAHIT on Biowulf';
MEGAHIT on Biowulf


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



MEGAHIT is a single node assembler for large and complex metagenomics NGS reads, such as soil. It makes use of succinct de Bruijn graph (SdBG) to achieve low memory assembly.



### References:


* Li, D., Luo, R., Liu, C.M., Leung, C.M., Ting, H.F., Sadakane, K., Yamashita, H. and Lam, T.W., 2016. **MEGAHIT v1.0: A Fast and Scalable Metagenome Assembler driven by Advanced Methodologies and Community Practices.** Methods. doi: [10.1016/j.ymeth.2016.02.020](https://doi.org/10.1016/j.ymeth.2016.02.020)
* Li, D., Liu, C-M., Luo, R., Sadakane, K., and Lam, T-W., (2015) **MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph.** Bioinformatics, doi: [10.1093/bioinformatics/btv033](https://doi.org/10.1093/bioinformatics/btv033)


Documentation
* [MEGAHIT Main Site](https://github.com/voutcn/megahit)
* [MEGAHIT wiki](https://github.com/voutcn/megahit/wiki)


Important Notes
* Module Name: megahit (see [the modules page](/apps/modules.html) for more information)
 * Be sure to set the max memory, threads, and [local scratch](/docs/userguide.html#local) temporary directory. See the example interactive job below.
* Environment variables set 
	+ MEGAHIT\_HOME* Example files in $MEGAHIT\_HOME/share/megahit/test\_data



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --mem 8g --cpus-per-task 2 --gres lscratch:10**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job
[user@cn3144 ~]$ **module load megahit**
[+] Loading megahit, version 1.1.4... 
[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144 46116226]$ **cp -a $MEGAHIT\_HOME/share/megahit/test\_data/\* .**
[user@cn3144 46116226]$ **SLURM\_MEM\_PER\_NODE\_BYTES="${SLURM\_MEM\_PER\_NODE}000000"**
[user@cn3144 46116226]$ **megahit --12 r1.il.fa.gz --memory ${SLURM\_MEM\_PER\_NODE\_BYTES} `# max memory to use` -t $SLURM\_CPUS\_PER\_TASK --tmp-dir /lscratch/$SLURM\_JOB\_ID**
--- [Fri Feb  1 19:59:29 2019] Start assembly. Number of CPU threads 2 ---
--- [Fri Feb  1 19:59:29 2019] Available memory: 270168522752, used: 8192000000
--- [Fri Feb  1 19:59:29 2019] Converting reads to binaries ---
    [read_lib_functions-inl.h  : 209]     Lib 0 (readsInterleaved1.fa.gz): interleaved, 200 reads, 100 max length
    [utils.h                   : 126]     Real: 0.0035	user: 0.0010	sys: 0.0029	maxrss: 6964
--- [Fri Feb  1 19:59:30 2019] k-max reset to: 119 ---
--- [Fri Feb  1 19:59:30 2019] k list: 21,29,39,59,79,99,119 ---
--- [Fri Feb  1 19:59:30 2019] Extracting solid (k+1)-mers for k = 21 ---
--- [Fri Feb  1 19:59:30 2019] Building graph for k = 21 ---
--- [Fri Feb  1 19:59:30 2019] Assembling contigs from SdBG for k = 21 ---
--- [Fri Feb  1 19:59:30 2019] Local assembling k = 21 ---
--- [Fri Feb  1 19:59:30 2019] Extracting iterative edges from k = 21 to 29 ---
--- [Fri Feb  1 19:59:30 2019] Building graph for k = 29 ---
--- [Fri Feb  1 19:59:30 2019] Assembling contigs from SdBG for k = 29 ---
--- [Fri Feb  1 19:59:30 2019] Local assembling k = 29 ---
--- [Fri Feb  1 19:59:30 2019] Extracting iterative edges from k = 29 to 39 ---
--- [Fri Feb  1 19:59:30 2019] Building graph for k = 39 ---
--- [Fri Feb  1 19:59:30 2019] Assembling contigs from SdBG for k = 39 ---
--- [Fri Feb  1 19:59:30 2019] Local assembling k = 39 ---
--- [Fri Feb  1 19:59:31 2019] Extracting iterative edges from k = 39 to 59 ---
--- [Fri Feb  1 19:59:31 2019] Building graph for k = 59 ---
--- [Fri Feb  1 19:59:31 2019] Assembling contigs from SdBG for k = 59 ---
--- [Fri Feb  1 19:59:31 2019] Local assembling k = 59 ---
--- [Fri Feb  1 19:59:31 2019] Extracting iterative edges from k = 59 to 79 ---
--- [Fri Feb  1 19:59:31 2019] Building graph for k = 79 ---
--- [Fri Feb  1 19:59:31 2019] Assembling contigs from SdBG for k = 79 ---
--- [Fri Feb  1 19:59:31 2019] Local assembling k = 79 ---
--- [Fri Feb  1 19:59:31 2019] Extracting iterative edges from k = 79 to 99 ---
--- [Fri Feb  1 19:59:31 2019] Merging to output final contigs ---
--- [STAT] 1 contigs, total 1207 bp, min 1207 bp, max 1207 bp, avg 1207 bp, N50 1207 bp
--- [Fri Feb  1 19:59:31 2019] ALL DONE. Time elapsed: 1.405510 seconds ---
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. megahit.sh). For example:



```

#!/bin/bash
set -e
module load megahit

test -n "$SLURM_CPUS_PER_TASK" || SLURM_CPUS_PER_TASK=2
test -n "$SLURM_MEM_PER_NODE" || SLURM_MEM_PER_NODE=1500

SLURM_MEM_PER_NODE_BYTES=${SLURM_MEM_PER_NODE}000000

megahit --12 readsInterleaved1.fa.gz --memory ${SLURM_MEM_PER_NODE_BYTES} -t $SLURM_CPUS_PER_TASK --tmp-dir /lscratch/$SLURM_JOB_ID

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] --gres lscratch:10 megahit.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. megahit.swarm). For example:



```

megahit --12 readsInterleaved1.fa.gz --memory ${SLURM_MEM_PER_NODE}000000 -t $SLURM_CPUS_PER_TASK --tmp-dir /lscratch/$SLURM_JOB_ID
megahit --12 readsInterleaved2.fa.gz --memory ${SLURM_MEM_PER_NODE}000000 -t $SLURM_CPUS_PER_TASK --tmp-dir /lscratch/$SLURM_JOB_ID
megahit --12 readsInterleaved3.fa.gz --memory ${SLURM_MEM_PER_NODE}000000 -t $SLURM_CPUS_PER_TASK --tmp-dir /lscratch/$SLURM_JOB_ID
megahit --12 readsInterleaved4.fa.gz --memory ${SLURM_MEM_PER_NODE}000000 -t $SLURM_CPUS_PER_TASK --tmp-dir /lscratch/$SLURM_JOB_ID

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f megahit.swarm -g # -t # --gres lscratch:# --module megahit
```

where


|  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --gres lscratch:# Number of gigabytes of [local scratch space](/docs/userguide.html#local) for each subjob in the swarm.
 | --module megahit Loads the MEGAHIT module for each subjob in the swarm
 | |
 | |
 | |
 | |








