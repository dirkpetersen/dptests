

document.querySelector('title').textContent = 'foldseek on Biowulf';
foldseek on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |




> 
>  Foldseek enables fast and sensitive comparisons of large structure
>  sets. It reaches sensitivities similar to state-of-the-art structural
>  aligners while being at least 20,000 times faster.
>  


Foldseek documentation

### References:


* M. van Kempen, S. Kim, C. Tumescheit, M. Mirdita, J. Söding, and M. Steinegger.
 *Foldseek: fast and accurate protein structure search.* BioRxiv, 2022
 [doi:10.1101/2022.02.07.479398](https://www.biorxiv.org/content/10.1101/2022.02.07.479398")


Documentation
* foldseek on[GitHub](https://github.com/steineggerlab/foldseek)
* foldseek on the web at [search.foldseek.com](https://search.foldseek.com/)


Important Notes
* Module Name: foldseek (see [the modules page](/apps/modules.html) for more information)
* foldseek is multithreaded and not slurm aware. Please specify threads to match the number of allocated CPUs
 to avoid overloaded jobs
* Example files in `$FOLDSEEK_TEST_DATA`
* Reference data in $FOLDSEEK\_DB/



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) for a short tutorial:



```

[user@biowulf]$ **sinteractive --mem=24g --cpus-per-task=12 --gres=lscratch:20**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load foldseek/5-53465f0**

```

Create a database from PDB or mmCif files to speed up repeated searches



```

[user@cn3144]$ **foldseek createdb --threads $SLURM\_CPUS\_PER\_TASK $FOLDSEEK\_TEST\_DATA testDB**
MMseqs Version:         3c64211f59830702e2a369d6e6f0d8a7492c27fa
Chain name mode         0
Write lookup file       1
Threads                 12
Verbosity               3

Output file: testDB
[=================================================================] 100.00% 26 0s 57ms
Time for merging to testDB_ss: 0h 0m 0s 0ms
Time for merging to testDB_h: 0h 0m 0s 0ms
Time for merging to testDB_ca: 0h 0m 0s 0ms
Time for merging to testDB: 0h 0m 0s 0ms
Ignore 0 out of 26.
Too short: 0, incorrect  0.
Time for processing: 0h 0m 0s 94ms

[user@cn3144]$ **foldseek easy-search --threads $SLURM\_CPUS\_PER\_TASK \
 $FOLDSEEK\_TEST\_DATA/d1asha\_ testDB aln.m8 /lscratch/$SLURM\_JOB\_ID/tmp**
...
[user@cn3144]$ **head aln.m8**
d1asha_ d1asha_ 1.000   147     0       0       1       147     1       147     2.859E-22       1061
d1asha_ d1x9fd_ 0.173   143     111     0       3       145     5       139     2.716E-05       265
d1asha_ d2w72b_ 0.182   145     112     0       1       145     4       141     4.463E-04       208
d1asha_ d1itha_ 0.131   145     119     0       1       145     3       140     6.611E-04       200
d1asha_ d1mbaa_ 0.118   143     120     0       4       146     6       142     6.611E-04       200

```

There are also pre-build foldseek databases for AlphafoldDB and PDB in $FOLDSEEK\_DB.



```

[user@cn3144]$ **foldseek easy-search --threads $SLURM\_CPUS\_PER\_TASK \
 --format-mode 4 \
 --format-output query,target,taxid,taxname,fident,alnlen,mismatch,qstart,qend,tstart,tend,evalue \
 $FOLDSEEK\_TEST\_DATA/d1asha\_ $FOLDSEEK\_DB/Alphafold-SwissProt aln\_afdb.tsv /lscratch/$SLURM\_JOB\_ID/tmp**
[user@cn3144]$

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```

In the following example individual steps are carried out separately and
hits are re-aligned with TM-align. As an alternative `--alignment-type
1` for easy-search uses TM-align and write the TMscore into the evalue
field. TM-score is in the range of (0,1]. 1 indicates a perfect match between
two structures. Scores below 0.17 correspond to randomly chosen unrelated
proteins whereas structures with a score higher than 0.5 assume generally the
same fold. The .tsv formatted outputs below can be merged to combine results.



```

[user@cn3144]$ **foldseek createdb --threads $SLURM\_CPUS\_PER\_TASK $FOLDSEEK\_TEST\_DATA queryDB**
[user@cn3144]$ **foldseek search -s 9.5 -a --threads $SLURM\_CPUS\_PER\_TASK queryDB $FOLDSEEK\_DB/PDB alnDB /lscratch/$SLURM\_JOB\_ID/tmp**
[user@cn3144]$ **foldseek convertalis --threads $SLURM\_CPUS\_PER\_TASK queryDB $FOLDSEEK\_DB/PDB alnDB aln.tsv**
[user@cn3144]$ **foldseek aln2tmscore --threads $SLURM\_CPUS\_PER\_TASK queryDB $FOLDSEEK\_DB/PDB alnDB alntmDB**
[user@cn3144]$ **foldseek createtsv --threads $SLURM\_CPUS\_PER\_TASK queryDB $FOLDSEEK\_DB/PDB alntmDB aln\_tm.tsv**
[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. foldseek.sh), which uses the input file 'foldseek.in'. For example:



```

#!/bin/bash
module load foldseek/1-2dd3b2f
foldseek easy-search --threads $SLURM_CPUS_PER_TASK \
    --format-mode 4 \
    --format-output query,target,taxid,taxname,fident,alnlen,mismatch,qstart,qend,tstart,tend,evalue
    $FOLDSEEK_TEST_DATA/d1asha_ $FOLDSEEK_DB/Alphafold-SwissProt aln_afdb.tsv /lscratch/$SLURM_JOB_ID/tmp

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=6 --mem=6G foldseek.sh
```







