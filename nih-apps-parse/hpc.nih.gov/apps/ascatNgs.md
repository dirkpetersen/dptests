

document.querySelector('title').textContent = 'Ascat NGS on Biowulf';
Ascat NGS on Biowulf


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



AscatNGS contains the Cancer Genome Projects workflow implementation of the ASCAT copy number algorithm for paired end sequencing.



### References:


* [Raine, Keiran M., et al. "ascatNgs: Identifying somatically acquired copy‐number alterations from whole‐genome sequencing data." *Current protocols in bioinformatics* 56.1 (2016): 15-9.](https://currentprotocols.onlinelibrary.wiley.com/doi/full/10.1002/cpbi.17)


Documentation
* [Ascat NGS Main Site](https://github.com/cancerit/ascatNgs)
* [Ascat NGS Wiki](https://github.com/cancerit/ascatNgs/wiki)


Important Notes
* Module Name: ascatngs (see [the modules page](/apps/modules.html) for more information)
* Ascat NGS is part of the Cancer Genome Project and is closely related to the programs [BRASS](/apps/BRASS.html) and [cgpBattenberg](/apps/cgpBattenberg.html) as well as the utilites [VAGrENT](https://github.com/cancerit/VAGrENT) and [PCAP-core](https://github.com/cancerit/PCAP-core). All of these programs can be added to your path using the cancerit-wgs module. To get the most recent versions of all of these, use the cancerit-wgs/latest module version. 
* Multithreaded app (use -c option)



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf ~]$ **sinteractive -c4 --mem=8g --gres=lscratch:10**
salloc.exe: Pending job allocation 11188180
salloc.exe: job 11188180 queued and waiting for resources
salloc.exe: job 11188180 has been allocated resources
salloc.exe: Granted job allocation 11188180
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0880 are ready for job
srun: error: x11: no local DISPLAY defined, skipping
error: unable to open file /tmp/slurm-spank-x11.11188180.0
slurmstepd: error: x11: unable to read DISPLAY value

[user@cn0880 ~]$ **module load ascatngs**
[+] Loading ascatngs  4.5.0  on cn0880
[+] Loading singularity  3.7.2  on cn0880

[user@cn0880 ~]$ **ascat.pl**

ERROR: Option must be defined.

Usage:
    ascat.pl [options]

      Please define as many of the parameters as possible

      Required parameters

        -outdir       -o    Folder to output result to.
        -tumour       -t    Tumour BAM/CRAM/counts file (counts must be .gz)
        -normal       -n    Normal BAM/CRAM/counts file (counts must be .gz)
        -reference    -r    Reference fasta
        -snp_gc       -sg   Snp GC correction file
        -protocol     -pr   Sequencing protocol (e.g. WGS, WXS)
        -gender       -g    Sample gender (XX, XY, L, FILE)
                              For XX/XY see '-gc'
                              When 'L' see '-l'
                              FILE - matched normal is_male.txt from ascatCounts.pl

      Targeted processing (further detail under OPTIONS):
        -process      -p    Only process this step then exit, optionally set -index
        -index        -i    Optionally restrict '-p' to single job
        -limit        -x    Specifying 2 will balance processing between '-i 1 & 2'
                            Must be paired with '-p allele_count'

      Optional parameters
        -genderChr    -gc   Specify the 'Male' sex chromosome: Y,chrY...
        -species      -rs   Reference species [BAM HEADER]
        -assembly     -ra   Reference assembly [BAM HEADER]
        -platform     -pl   Seqeuncing platform [BAM HEADER]
        -minbasequal  -q    Minimum base quality required before allele is used. [20]
        -cpus         -c    Number of cores to use. [1]
                            - recommend max 2 during 'input' process.
        -locus        -l    Using a list of loci, default when '-L' [share/gender/GRCh37d5_Y.loci]
                            - these are loci that will not be present at all in a female sample
        -force        -f    Force completion - solution not possible
                            - adding this will result in successful completion of analysis even
                              when ASCAT can't generate a solution.  A default copynumber of 5/2
                              (tumour/normal) and contamination of 30% will be set along with a
                              comment in '*.samplestatistics.csv' to indicate this has occurred.
        -purity       -pu   Purity (rho) setting for manual setting of sunrise plot location
        -ploidy       -pi   Ploidy (psi) setting for manual setting of sunrise plot location
        -noclean      -nc   Finalise results but don't clean up the tmp directory.
                            - Useful when including a manual check and restarting ascat with new pu and pi params.
        -nobigwig     -nb   Don't generate BigWig files.
        -t_name       -tn   Tumour name to use when using count files as input
        -n_name       -nn   Noraml name to use when using count files as input

      Other
        -help         -h    Brief help message
        -man          -m    Full documentation.
        -version      -v    Ascat version number


[user@cn0880 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn0880 11188180]$ **cp /path/to/tumor.bam .**

[user@cn0880 11188180]$ **cp /path/to/normal.bam .**

[user@cn0880 11188180]$ **ascat.pl -o output -t tumor.bam -n normal.bam \
 -r /fdb/igenomes/Homo\_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa \
 -snp\_gc SnpGcCorrections.tsv -pr wgs -g XX -c $SLURM\_CPUS\_PER\_TASK**
[snip...]

[user@cn0880 11188180]$ **exit**
exit
srun: error: cn0880: task 0: Exited with exit code 2
salloc.exe: Relinquishing job allocation 11188180

[user@biowulf ~]$

```


Note: You need to generate the SnpGcCorrections.tsv (as can be found in: https://github.com/cancerit/ascatNgs/wiki/Convert-SnpPositions.tsv-to-SnpGcCorrections.tsv) or downloaded (https://github.com/cancerit/ascatNgs/wiki/Human-reference-files-from-1000-genomes-VCFs).

Generates LogR.txt and BAF.txt, which can be used to generate non-segmented plots (see https://www.crick.ac.uk/peter-van-loo/software/ASCAT)




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. AscatNGS.sh). For example:



```

#!/bin/bash
module load ascatngs
export GENOME=/fdb/igenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/
ascat.pl -o output -t tumor.bam -n normal.bam -r $GENOME/genome.fa -snp_gc SnpGcCorrections.tsv -pr wgs -g XX -c $SLURM_CPUS_PER_TASK

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=16 --mem=30g AscatNGS.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. AscatNGS.swarm). For example:



```

ascat.pl -o output1 -t tumor1.bam -n normal1.bam ... -c $SLURM_CPUS_PER_TASK
ascat.pl -o output2 -t tumor2.bam -n normal2.bam ... -c $SLURM_CPUS_PER_TASK
ascat.pl -o output3 -t tumor3.bam -n normal3.bam ... -c $SLURM_CPUS_PER_TASK

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f AscatNGS.swarm -g 30 -t 16 --module ascatngs
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module ascatngs Loads the Ascat NGS module for each subjob in the swarm 
 | |
 | |
 | |








