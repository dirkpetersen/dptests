

document.querySelector('title').textContent = 'pplacer on Biowulf';
pplacer on Biowulf


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



Pplacer places query sequences on a fixed reference phylogenetic tree to maximize phylogenetic likelihood or posterior probability 
according to a reference alignment. Pplacer is designed to be fast, to give useful information about uncertainty, and to offer advanced 
visualization and downstream analysis.



### References:


* pplacer was developed by the Matsen group at the Fred Hutchinson Cancer Research Center. See the [pplacer website](http://matsen.fhcrc.org/pplacer/) 
for references.


Documentation
* [Pplacer documentation](http://matsen.fhcrc.org/pplacer/) at the Matsen lab, FHCC


Important Notes
* Module Name: pplacer (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded
* Example files in /usr/local/apps/pplacer/tutorial.tar.gz



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

[user@cn3144 ~]$ **mkdir /data/$USER/pplacer; cd /data/$USER/pplacer**

[user@cn3144 pplacer]$ **tar xvzf /usr/local/apps/pplacer/tutorial.tar.gz**
[...]

[user@cn3144 pplacer]$ **cd fhcrc-microbiome-demo-730d268**

[user@cn3144 pplacer]$ **module avail pplacer**

----------------------------------------- /usr/local/lmod/modulefiles --------------------------
   pplacer/1.1

[user@cn3144 pplacer]$ **module load pplacer**
[+] Loading pplacer 1.1

[user@cn3144 pplacer]$ **sh pplacer\_demo.sh**
# Phylogenetic placement
# ----------------------

# This makes p4z1r2.jplace, which is a "place" file in JSON format.  Place files
# contain information about collections of phylogenetic placements on a tree.
# You may notice that one of the arguments to this command is
# `vaginal_16s.refpkg`, which is a "reference package". Reference packages are
# simply an organized collection of files including a reference tree, reference
# alignment, and taxonomic information. We have the beginnings of a
# [database](http://microbiome.fhcrc.org/apps/refpkg/) of reference packages
# and some [software](http://github.com/fhcrc/taxtastic) for putting them
# together.
pplacer -c vaginal_16s.refpkg src/p4z1r36.fasta
Running pplacer v1.1.alpha19-0-g807f6f3 analysis on src/p4z1r36.fasta...
Found reference sequences in given alignment file. Using those for reference alignment.
Pre-masking sequences... sequence length cut from 2196 to 275.
Determining figs... figs disabled.
Allocating memory for internal nodes... done.
Caching likelihood information on reference tree... done.
Pulling exponents... done.
Preparing the edges for baseball... done.
warning: rank below_species not represented in the lineage of any sequence in reference package vaginal_16s.
[...]
pause
Please press return to continue...

![](/images/pplacer.jpg)

# That's it for the demo. For further information, please consult the
# [pplacer documentation](http://matsen.github.com/pplacer/).
echo "Thanks!"
Thanks!

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. pplacer.sh). For example:



```

#!/bin/bash
set -e
module load pplacer
mkdir /data/$USER/pplacer
cd /data/$USER/pplacer
tar xvzf /usr/local/apps/pplacer/tutorial.tar.gz
cd fhcrc-microbiome-demo-730d268

pplacer -c vaginal_16s.refpkg src/p4z1r36.fasta

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch  [--mem=#] pplacer.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. pplacer.swarm). For example:



```

guppy kr_heat -c vaginal_16s.refpkg/ file1a.jplace file2a.jplace
guppy kr_heat -c vaginal_16s.refpkg/ file1b.jplace file2b.jplace
guppy kr_heat -c vaginal_16s.refpkg/ file1c.jplace file2c.jplace

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f pplacer.swarm [-g #]  --module pplacer
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --module pplacer Loads the pplacer module for each subjob in the swarm 
 | |
 | |








