

document.querySelector('title').textContent = 'Scramble: a tool for mobile element insertion detection';
**Scramble: a tool for mobile element insertion detection**


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



Scramble is a mobile element insertion (MEI) detection tool. It identifies clusters of soft clipped reads
in a BAM file, builds consensus sequences, aligns to representative L1Ta, AluYa5, and SVA-E sequences, and
outputs MEI calls. 



### References:


* Rebecca I. Torene, Kevin Galens, Shuxi Liu, Kevin Arvai, Carlos Borroto, Julie Scuffins, 
Zhancheng Zhang, Bethany Friedman, Hana Sroka, Jennifer Heeley, Erin Beaver, Lorne Clarke, 
Sarah Neil, Jagdeep Walia, Danna Hull, Jane Juusola, and Kyle Retterer.  

*Mobile element insertion detection in 89,874 clinical exomes*   

[Genetics in Medicine](https://www.nature.com/articles/s41436-020-0749-x) (2020).


Documentation
* [Scramble GitHub page](https://github.com/GeneDx/scramble)


Important Notes
* Module Name: Scramble (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Environment variables set
	+ **SCRAMBLE\_HOME**  Scramble installation directory
	+ **SCRAMBLE\_BIN**      Scramble executable directory
	+ **SCRAMBLE\_DATA**  Scramble sample data directory
	+ If you are using your own reference file make sure you generate **\*.nhr**, **\*.nin**, and **\*.nsq** files using **makeblastd** as follows:  
	
	*module load load ncbi-toolkit  
	
	 makeblastdb -in file.fasta -input\_type fasta -dbtype nucl*



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf ~]$ **sinteractive --mem=4g**
salloc.exe: Pending job allocation 56730292
salloc.exe: job 56730292 queued and waiting for resources
salloc.exe: job 56730292 has been allocated resources
salloc.exe: Granted job allocation 56730292
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3148 are ready for job
[user@cn3148 ~]$ **module load Scramble** 
[+] Loading singularity  3.5.3  on cn3148
[+] Loading Scramble 0.0.20190211.82c78b9  ...

```

Copy sample data into your current directory:

```

[user@cn3148 ~]$ **cp $SCRAMBLE\_DATA/\* .**

```
 
You can run Scramble in two different ways.   


