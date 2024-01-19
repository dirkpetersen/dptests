

document.querySelector('title').textContent = ' weblogo on Biowulf2 &amp; Helix';

weblogo on Biowulf2 & Helix



|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Interactive job](#int)
[Batch job](#serial)
[Swarm of jobs](#swarm)
 |



The weblogo package includes a tool ( `seqlogo`) to generate
sequence logos from aligned sequences in fastq format.



Note that weblogo changed significantly from version 2.X to 3.X. The commands
listed here are for versions 3.X


### References


* G. E. Crooks, G. Hon, J. M. Chandonia, and S. E. Brenner. *WebLogo: A sequence 
 logo generator.* Genome Research 2004, 14:1188-1190. 
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/15173120) | 
 [PMC](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC419797/) | 
 [Journal](http://genome.cshlp.org/content/14/6/1188.long)


### Documentation


* [Weblogo2](http://weblogo.berkeley.edu/)


Important Notes
* Module Name: weblogo (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded
* Example files in $WEBLOGO\_DATA



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cann
ot be run as batch jobs.
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

[user@cn3144 ~]$ **module load weblogo**
[+] Loading weblogo 3.6
[user@cn3144 ~]$ **cp $WEBLOGO\_DATA/\* .**
[user@cn3144 ~]$ **weblogo -f ctcf.fa -F png -o ctcf.png**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

Visualization of the output image ctcf.png:
![](/images/weblogo_fig1.png)


Batch job on Biowulf2
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch script similar to the following example:



```

#! /bin/bash
#SBATCH --mem=100m
#SBATCH --time=5
# this file is seqlogo_job.sh

module load weblogo || exit 1
cp $WEBLOGO_DATA/* .
weblogo -f ctcf.fa -F png -o ctcf.png


```

Submit to the queue with [sbatch](/docs/userguide.html):



```

biowulf$ **sbatch seqlogo\_job.sh**

```



Swarm of jobs on Biowulf2
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical
 resources.
Create a swarm command file similar to the following example:



```

# this file is seqlogo_jobs.swarm
seqlogo -f site1.fa -F eps -o site1 -abcMnY
seqlogo -f site2.fa -F eps -o site2 -abcMnY
seqlogo -f site3.fa -F eps -o site3 -abcMnY
seqlogo -f site4.fa -F eps -o site4 -abcMnY

```

And submit to the queue with [swarm](/apps/swarm.html)



```

biowulf$ **swarm -f seqlogo\_jobs.swarm --time=5**

```







