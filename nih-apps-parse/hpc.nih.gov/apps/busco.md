

document.querySelector('title').textContent = 'busco on Biowulf';
busco on Biowulf


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



BUSCO completeness assessments employ sets of Benchmarking Universal Single-Copy Orthologs from [OrthoDB](http://www.orthodb.org) to provide quantitative measures of the completeness of genome assemblies, annotated gene sets, and transcriptomes in terms of expected gene content. 



### References:


* Robert M. Waterhouse, Mathieu Seppey, Felipe A. Simao, Mose Manni, Panagiotis Ioannidis, Guennadi Klioutchnikov, Evgenia V. Kriventseva, and Evgeny M. Zdobnov. [BUSCO applications from quality assessments to gene prediction and phylogenomics.](https://doi.org/10.1093/molbev/msx319)
*Mol Biol Evol, published online Dec 6, 2017*
* Felipe A. Simao, Robert M. Waterhouse, Panagiotis Ioannidis, Evgenia V. Kriventseva, and Evgeny M. Zdobnov. [BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs.](https://doi.org/10.1093/bioinformatics/btv351)
*Bioinformatics, published online June 9, 2015*


 **BUSCO configuration files between different BUSCO versions will likely be incompatible.**  

 Please rename the configuration file in your home directory, busco.config to, say, busco.config.OLD. When you load different versions of BUSCO, a new configuration file will be created in your home directory ~/busco.config but ***only if*** it doesn't already exist. As a result it is important to rename previous configuration files when switching between different BUSCO versions. You can also set the contents of the environment variable BUSCO\_CONFIG\_FILE to a path and file name of your choosing if you are planning to go back and forth between different versions of BUSCO.



Documentation
* [BUSCO Main Site](https://gitlab.com/ezlab/busco)* [BUSCO User Guide](https://busco.ezlab.org/busco_userguide.html)


Important Notes
* Module Name: busco (see [the modules page](/apps/modules.html) for more information)
* BUSCO versions 4 and greater perform automated download of all necessary files and datasets to conduct a run.
* BUSCO requires a configuration file that sets the paths for required executables (tblastn, makeblastdb, augustus, hmmsearch) as well as parameters that the user wants to set for their busco runs. When you load the busco module, the default configuration file is copied to /home/$USER/busco.config. The paths for tblastn, makeblastdb, augustus, hmmsearch are already set correctly in this file. You can edit the file and set any additional busco parameters that you want. 
* [Augustus](https://github.com/Gaius-Augustus/Augustus) is a gene prediction program for eukaryotes which is required by BUSCO. Augustus requires a writable configuration directory. When you load the busco module, the Augustus config directory will be copied to /home/$USER/augustus.config, and Augustus runs will update this directory. 
* Environment variables set 
	+ BUSCO\_CONFIG\_FILE - set to /home/$USER/busco.config
	+ AUGUSTUS\_CONFIG\_PATH - set to /home/$USER/augustus.config
	+ BUSCO\_EXAMPLE\_DATA - set to /usr/local/apps/busco/sample\_data



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive -c 6 --mem=6G**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load busco**
$ module load busco
[+] Loading singularity  3.7.1  on cn3144
[+] Loading busco  4.1.3
Busco Config file /home/user/busco.config  does not exist
   copying from /usr/local/apps/busco/4.1.3/share/config.ini.v4.1.3
 You should edit this file to set your own parameters
Augustus Config directory /home/user/augustus.config  does not exist
   copying from /usr/local/apps/busco/4.1.3/share/augustus_config

*#####Edit the busco.config file to set any additional parameters desired*

[user@cn3144 ~]$ **nano ~/busco.config**

*##### Copy the sample data set* 

[user@cn3144 ~]$ **cp $BUSCO\_EXAMPLE\_DATA/bacteria/genome.fna .**

*##### Run BUSCO*
[user@cn3144 ~]$ **busco -i genome.fna -c 6 -m geno -f --out test\_bacteria**
INFO:   ***** Start a BUSCO v4.1.3 analysis, current time: 03/09/2021 12:22:58 *****
INFO:   Configuring BUSCO with /home/user/busco.config
INFO:   Mode is genome
INFO:   Input file is genome.fna
INFO:   Downloading information on latest versions of BUSCO data...
WARNING:        Running Auto Lineage Selector as no lineage dataset was specified. This will take a little longer than normal. If you know what lineage dataset you want to use, please specify this in the config file or using the -l (--lineage-dataset) flag in the command line.
INFO:   No lineage specified. Running lineage auto selector.

INFO:   ***** Starting Auto Select Lineage *****
        This process runs BUSCO on the generic lineage datasets for the domains archaea, bacteria and eukaryota. Once the optimal domain is selected, BUSCO automatically attempts to find the most appropriate BUSCO dataset to use based on phylogenetic placement.
        --auto-lineage-euk and --auto-lineage-prok are also available if you know your input assembly is, or is not, an eukaryote. See the user guide for more information.
        A reminder: Busco evaluations are valid when an appropriate dataset is used, i.e., the dataset belongs to the lineage of the species to test. Because of overlapping markers/spurious matches among domains, busco matches in another domain do not necessarily mean that your genome/proteome contains sequences from this domain. However, a high busco score in multiple domains might help you identify possible contaminations.
INFO:   Downloading file 'https://busco-data.ezlab.org/v4/data/lineages/archaea_odb10.2020-03-06.tar.gz'
INFO:   Decompressing file '/home/user/busco_downloads/lineages/archaea_odb10.tar.gz'
INFO:   Running BUSCO using lineage dataset archaea_odb10 (prokaryota, 2020-03-06)
INFO:   ***** Run Prodigal on input to predict and extract genes *****
INFO:   Running Prodigal with genetic code 11 in single mode
INFO:   Running 1 job(s) on prodigal, starting at 03/09/2021 12:23:07
INFO:   [prodigal]      1 of 1 task(s) completed
INFO:   Genetic code 11 selected as optimal
INFO:   ***** Run HMMER on gene sequences *****
[...]


[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. busco.sh). The example below assumes that you have downloaded the lineage-specific profile library and edited the BUSCO config as desired. 


```

#!/bin/bash
set -e
module load busco/4.1.3
busco -m genome -i target.fa -o test2 -l bacteria_odb10

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] busco.sh
```









