

document.querySelector('title').textContent = 'mash on Biowulf';
mash on Biowulf


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



Mash uses MinHash hashing to reduce large sequences to a representative sketch.
Distances between sketches of sequences can be calculated very rapidly and
can provide an estimate of average nucleotide identity. Sketches of all the
genomes in RefSeq 70 are only ~90MB (at a kmer size of 16 using 400 hashes).



### References:


* Brian D. Ondov *et al.* *Fast genome 
 and metagenome distance estimation using MinHash*. bioRxiv 2015:
 http://dx.doi.org/10.1101/029827
 [bioRxiv](http://biorxiv.org/content/early/2015/10/26/029827.article-info)


Documentation
* [GitHub](https://github.com/marbl/Mash)
* [Docs](http://mash.readthedocs.org/en/latest/index.html)


Important Notes
* Module Name: mash (see [the modules page](/apps/modules.html) for more information)
* Example files in can be found in `$MASH_TEST_DATA`
* Some mash tools can be run with multiple threads.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --gres=lscratch:20**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load mash**
[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

```

Some example data is included in the `TEST_DATA` directory
inside the application directory.



```

[user@cn3144 ~]$ **ls $MASH\_TEST\_DATA**
BA000007.2.fna  CP009273.1.fna  NC_000913.2.fna  RefSeqGenomes_V70.msh  SRR292770_1.fastq.gz
[user@cn3144 ~]$ **cp $MASH\_TEST\_DATA/\* .**

```

Estimate the distance between two *E. coli* genomes



```

[user@cn3144 ~]$ **mash dist BA000007.2.fna NC\_000913.2.fna**
/usr/local/apps/mash/TEST_DATA/BA000007.2.fna   /usr/local/apps/mash/TEST_DATA/NC_000913.2.fna
0.0222766      0       456/1000

```

The result shows the reference sequence, the query sequence, the distance
estimate, the p value, and the number of matching hashes.


Instead of calculating sketches each time they can be precalculated. For example,
we can sketch the two genomes from above



```

[user@cn3144 ~]$ **mash sketch -o 2coli.msh BA000007.2.fna NC\_000913.2.fna**
Sketching /usr/local/apps/mash/TEST_DATA/BA000007.2.fna...
Sketching /usr/local/apps/mash/TEST_DATA/NC_000913.2.fna...
Writing to 2coli.msh...
[user@cn3144 ~]$ **ls -lh 2coli.msh**
-rw-r--r-- 1 user group  16K Dec  4 09:56 2coli.msh
[user@cn3144 ~]$ **mash info 2coli.msh**
Header:
  Hash function (seed):          MurmurHash3_x64_128 (42)
  K-mer size:                    21 (64-bit hashes)
  Alphabet:                      ACGT (canonical)
  Target min-hashes per sketch:  1000
  Sketches:                      2

Sketches:
  [Hashes]  [Length]  [ID]             [Comment]
  1000      5498450   BA000007.2.fna   -
  1000      4639675   NC_000913.2.fna  -


```

So the sketch for the 10MB genomes takes up 16kB. Now we can compare our 
query against the sketches with



```

[user@cn3144 ~]$ **mash dist 2coli.msh CP009273.1.fna**
BA000007.2.fna   CP009273.1.fna   0.0222766      0       456/1000
NC_000913.2.fna  CP009273.1.fna   0              0      1000/1000

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. mash.sh), which uses the input file 'mash.in'. For example:



```

#! /bin/bash

function fail {
  echo "$@" >&2
  exit 1
}

module load mash/2.0 || fail "could not load mash module"
mash sketch -i -k 21 -s 1000 -o coli.msh -l list_of_ecoli_genomes

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] mash.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. mash.swarm). For example:



```

mash sketch -o NC_000913.2.msh /usr/local/apps/mash/TEST_DATA/NC_000913.2.fna
mash sketch -o BA000007.2.msh /usr/local/apps/mash/TEST_DATA/BA000007.2.fna
mash sketch -o CP009273.1.msh /usr/local/apps/mash/TEST_DATA/CP009273.1.fna

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f mash.swarm [-g #] [-t #] --module mash
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module mash  Loads the mash module for each subjob in the swarm 
 | |
 | |
 | |








