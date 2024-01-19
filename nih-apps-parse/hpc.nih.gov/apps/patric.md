

document.querySelector('title').textContent = ' PATRIC';
 PATRIC


|  |
| --- |
| 
Quick Links
[Important notes](#notes)
[Download errors](#errors)
[Downloading data from PATRIC](#patric)
[Batch jobs](#sbatch)
[Swarm jobs](#swarm)
 |


The [PATRIC CLI](https://docs.patricbrc.org/cli_tutorial/index.html) allows programmatic access to the Pathosystems Resource Integration Center online database. 



Documentation
* [PATRIC BRC Web Database](https://patricbrc.org/)* [PATRIC CLI Tutorial](https://docs.patricbrc.org/cli_tutorial/index.html)
* [CLI Command Reference](https://docs.patricbrc.org/cli_tutorial/command_list/index.html)
* **Command line help:** Type the command followed by '-h'


Important Notes
* Module Name: patric (see [the modules page](/apps/modules.html) for more information)
* PATRIC commands are prefixed with p3-




Downloading data from PATRIC
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.

The primary use of PATRIC CLI tools is downloading and querying data from the PATRIC database. For instance, to query and download a specific baterial genome.:

```

[USER@biowulf]$ **sinteractive --gres=lscratch:40**
salloc: Pending job allocation 43361600
salloc: job 43361600 queued and waiting for resources
salloc: job 43361600 has been allocated resources
salloc: Granted job allocation 43361600
salloc: Waiting for resource configuration
salloc: Nodes cn2296 are ready for job

[USER@cn2296]$ **mkdir /data/$USER/patric\_genomes**

[USER@cn2296]$ **cd !$**

[USER@cn2296]$ **module load patric**

[USER@cn2296]$ **p3-all-genomes --eq genus,Paenibacillus -a genome\_name | head**
genome.genome_id        genome.genome_name
1333861.3       Paenibacillus forsythiae T98
1619309.3       Paenibacillus herberti strain R33
697284.3        Paenibacillus larvae subsp. larvae DSM 25430
697284.9        Paenibacillus larvae subsp. larvae DSM 25430 strain DSM 25430; ERIC_II
697286.3        Paenibacillus larvae subsp. larvae DSM 25719
147375.10       Paenibacillus larvae subsp. larvae strain Eric_V
1330551.3       Paenibacillus lemnae strain L7-75
1619311.3       Paenibacillus physcomitrellae strain XB
1414587.3       Paenibacillus polymyxa A18

[USER@cn2296]$ **p3-genome-fasta 1414587.3 > Paenibacillus\_polymyxa\_A18.fasta**

[USER@cn2296]$ **grep ">" Paenibacillus\_polymyxa\_A18.fasta**
>JWJJ01000003 scaffold
>JWJJ01000001 scaffold
>JWJJ01000005 scaffold
>JWJJ01000004 scaffold
>JWJJ01000022 scaffold
>JWJJ01000006 scaffold
>JWJJ01000026 scaffold
>JWJJ01000021 scaffold
>JWJJ01000024 scaffold
>JWJJ01000019 scaffold
>JWJJ01000014 scaffold
>JWJJ01000025 scaffold
>JWJJ01000012 scaffold
>JWJJ01000018 scaffold
>JWJJ01000016 scaffold
>JWJJ01000023 scaffold
>JWJJ01000017 scaffold
>JWJJ01000013 scaffold
>JWJJ01000011 scaffold
>JWJJ01000027 scaffold
>JWJJ01000008 scaffold
>JWJJ01000010 scaffold
>JWJJ01000002 scaffold
>JWJJ01000009 scaffold
>JWJJ01000007 scaffold
>JWJJ01000015 scaffold
>JWJJ01000029 scaffold
>JWJJ01000020 scaffold
>JWJJ01000028 scaffold

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. patric.sh). For example:



```

#!/bin/bash
set -e
module load patric
p3-genome-fasta 1414587.3 > 1414587.3.fasta

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=1 --mem=2g patric.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. patric.swarm). For example:



```

cd fasta_dir; p3-genome-fasta 1414587.3 > 1414587.3.fa
cd fasta_dir; p3-genome-fasta 697284.3 > 697284.3
cd fasta_dir; p3-genome-fasta 697286.3 > 697286.3
cd fasta_dir; p3-genome-fasta 1333861.3 > 1333861.3

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f patric.swarm [-g #] [-t #] --module patric
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module patric Loads the patric module for each subjob in the swarm
 | |
 | |
 | |












