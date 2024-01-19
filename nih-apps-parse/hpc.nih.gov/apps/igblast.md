

document.querySelector('title').textContent = 'igblast on Biowulf';
igblast on Biowulf


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



igblast is a tool developed at [NCBI](https://www.ncbi.nlm.nih.gov) for analysis of immunoglobulin variable domain sequences and T cell receptor (TR) sequences. In addition to performing a regular Blast search, igblast
* Reports the germline V, D and J gene matches to the query sequence.
 Annotates the immunoglobulin domains (FR1 through CDR3).
 Reveals the V(D)J junction details such as nucleotide homology between the ends of V(D)J segments and N nucleotide insertions.
 Indicates the whether the rearrangement is in-frame or out-of-frame.





Documentation
* [igblast docs](https://www.ncbi.nlm.nih.gov/igblast/intro.html)


Important Notes
* Module Name: igblast (see [the modules page](/apps/modules.html) for more information)
* Multithreaded
* Reference data in /fdb/igblast/



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

[user@cn3144 ~]$ **module load igblast samtools**

# Get the chromosome containing mouse Tcra (T-cell receptor alpha chain)
[user@cn3144 ~]$ **ln -s /fdb/igenomes/Mus\_musculus/NCBI/GRCm38/Sequence/Chromosomes/14.fa**

[user@cn3144 ~]$ **igblastn \
 -germline\_db\_V mouse\_gl\_V \
 -germline\_db\_J mouse\_gl\_J \
 -germline\_db\_D mouse\_gl\_D \
 -organism mouse \
 -domain\_system kabat \
 -query <(samtools faidx 14.fa 14:52427967-54224198) \
 -auxiliary\_data optional\_file/mouse\_gl.aux \
 -show\_translation \
 -outfmt 7** 

# IGBLASTN 2.6.1+
# Query: 14:52427967-54224198
# Database: mouse_gl_V mouse_gl_D mouse_gl_J
# Domain classification requested: kabat

# V-(D)-J rearrangement summary for query sequence (Top V gene match, Top J gene match, Chain type, stop codon, V-J frame, Productive, Strand).  Multiple equivalent top matches, if present, are separated by a comma.
		    21-11N/AVLYesN/ANo+

# V-(D)-J junction details based on top germline gene matches (V end, V-J junction, J start).  Note that possible overlapping nucleotides at VDJ junction (i.e, nucleotides that could be assigned to either rearranging gene) are indicated in parentheses (i.e., (TACT)) but are not included under the V, D, or J gene itself
		    CAAGGN/AN/A

# Hit table (the first field indicates the chain type of the hit)
# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, gaps, q. start, q. end, s. start, s. end, evalue, bit score
# 3 hits found
		    V14:52427967-5422419821-1189.15783524121582612159042853671.47e-1690.6
		    V14:52427967-5422419821-1186.905847245232845233642853673.74e-1586.0
		    V14:52427967-5422419821-1186.905847248342948343742853673.74e-1586.0

Total queries = 1
Total identifiable CDR3 = 0
Total unique clonotypes = 0

# BLAST processed 1 queries
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. igblast.sh). For example:



```

#!/bin/bash
set -e
 
cd /data/$USER/mydir
module load igblast samtools
ln -s /fdb/igenomes/Mus_musculus/NCBI/GRCm38/Sequence/Chromosomes/14.fa
igblastn \
  -germline_db_V mouse_gl_V \
  -germline_db_J mouse_gl_J \
  -germline_db_D mouse_gl_D \
  -organism mouse \
  -domain_system kabat \
  -query <(samtools faidx 14.fa 14:52427967-54224198) \
  -auxiliary_data optional_file/mouse_gl.aux \
  -show_translation \
  -outfmt 7
  -num_threads $SLURM_CPUS_PER_TASK

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command, e.g.



```
sbatch --cpus-per-task=8 --mem=10g igblast.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. igblast.swarm). For example:



```

igblastn   -germline_db_V mouse_gl_V   -germline_db_J mouse_gl_J   -germline_db_D mouse_gl_D \
  -organism mouse   -domain_system kabat   -query <(samtools faidx 14.fa 14:52427967-54224198) \
  -auxiliary_data optional_file/mouse_gl.aux   -show_translation   -outfmt 7
  -num_threads $SLURM_CPUS_PER_TASK
[...]

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f igblast.swarm [-g #] [-t #] --module igblast
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module igblast Loads the igblast module for each subjob in the swarm 
 | |
 | |
 | |








