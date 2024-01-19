

document.querySelector('title').textContent = 'Nirvana: clinical-grade annotation of genomic variants';
Nirvana: clinical-grade annotation of genomic variants


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



Nirvana provides clinical-grade annotation of genomic variants (SNVs, MNVs, insertions, deletions, indels, and SVs (including CNVs). It can be run as a stand-alone package or integrated into larger software tools that require variant annotation.



### References:


* Michael Stromberg, Rajat Roy, Julien Lajugie, Yu Jiang, Haochen Li, Elliott Margulies   

*Nirvana: Clinical Grade Variant Annotator*    

[Proc. of the 8th ACM Intern. Conf. on Bioinformatics, Computational Biology,and Health Informatics, August](https://dl.acm.org/doi/abs/10.1145/3107411.3108204) 2017 p.596, doi: https://doi.org/10.1145/3107411.3108204.


Documentation
* [Nirvana Github page](https://github.com/Illumina/Nirvana)
* [Nirvana Documentation](https://illumina.github.io/NirvanaDocumentation/)


Important Notes
* Module Name: nirvana (see [the modules page](/apps/modules.html) for more information)
* singlethreaded app
* Example files in /usr/local/apps/nirvana/TEST\_DATA* Reference data in /fdb/nirvana/



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=8g**
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load nirvana** 
[+] Loading nirvana  2.16.1 on cn3144
[user@cn3144 ~]$ **nirvana -h**
---------------------------------------------------------------------------
Nirvana                                             (c) 2021 Illumina, Inc.
Stromberg, Roy, Lajugie, Jiang, Li, and Kang                         3.16.1
---------------------------------------------------------------------------

USAGE: dotnet Nirvana.dll -i <vcf path> -c <cache prefix> --sd <sa dir> -r <ref path> -o <base output filenam>e
Annotates a set of variants

OPTIONS:
      --cache, -c <prefix>   input cache prefix
      --in, -i <t;path>        input VCF path
      --out, -o <file path>  output file path
      --ref, -r <path>       input compressed reference sequence path
      --sd <directory>       input supplementary annotation directory
      --force-mt             forces to annotate mitochondrial variants
      --disable-recomposition
                             don't recompose function relevant variants
      --legacy-vids          enables support for legacy VIDs
      --enable-dq            report DQ from VCF samples field
      --str <VALUE>          user provided STR annotation TSV file
      --help, -h             displays the help menu
      --version, -v          displays the version 

```

Download a sample VCF file:

```

[user@cn3144 ~]$  **curl -O https://illumina.github.io/NirvanaDocumentation/files/HiSeq.10000.vcf.gz** 

```

Run Nirvana on the sample file:

```

[user@cn3144 ~]$ **nirvana \
 -c $NIRVANA\_DATA/Cache/GRCh37/Both \
 --sd $NIRVANA\_DATA/SupplementaryAnnotation/GRCh37 \
 -r $NIRVANA\_DATA/References/Homo\_sapiens.GRCh37.Nirvana.dat \
 -i HiSeq.10000.vcf.gz \
 -o HiSeq.10000** 
---------------------------------------------------------------------------
Nirvana                                             (c) 2021 Illumina, Inc.
Stromberg, Roy, Lajugie, Jiang, Li, and Kang                         3.16.1
---------------------------------------------------------------------------

Initialization                                         Time     Positions/s
---------------------------------------------------------------------------
Cache                                               00:00:02.2
SA Position Scan                                    00:00:00.1       61,441

Reference                                Preload    Annotation   Variants/s
---------------------------------------------------------------------------
chr1                                    00:00:00.7  00:00:03.4        2,904

Summary                                                Time         Percent
---------------------------------------------------------------------------
Initialization                                      00:00:02.4       21.4 %
Preload                                             00:00:00.7        6.6 %
Annotation                                          00:00:03.4       30.2 %

Peak memory usage: 1.322 GB
Time: 00:00:10.8

[user@cn3144 ~]$ **nirvana \
 -c $NIRVANA\_DATA/Cache/GRCh38/Both \
 --sd $NIRVANA\_DATA/SupplementaryAnnotation/GRCh38 \
 -r $NIRVANA\_DATA/References/Homo\_sapiens.GRCh38.Nirvana.dat \
 -i HiSeq.10000.vcf.gz \
 -o HiSeq.10000** 
---------------------------------------------------------------------------
Nirvana                                             (c) 2021 Illumina, Inc.
Stromberg, Roy, Lajugie, Jiang, Li, and Kang                         3.16.1
---------------------------------------------------------------------------

Initialization                                         Time     Positions/s
---------------------------------------------------------------------------
Cache                                               00:00:02.1
SA Position Scan                                    00:00:00.1       66,470

Reference                                Preload    Annotation   Variants/s
---------------------------------------------------------------------------
chr1                                    00:00:01.3  00:00:02.9        3,325

Summary                                                Time         Percent
---------------------------------------------------------------------------
Initialization                                      00:00:02.3       19.0 %
Preload                                             00:00:01.3       11.3 %
Annotation                                          00:00:02.9       24.6 %

Peak memory usage: 1.346 GB
Time: 00:00:11.4
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. nirvana.sh). For example:



```

#!/bin/bash
set -e
module load nirvana
nirvana -c $NIRVANA_DATA/Cache/GRCh37/Both \
		--sd $NIRVANA_DATA/SupplementaryAnnotation/GRCh37 \
		-r $NIRVANA_DATA/References/Homo_sapiens.GRCh37.Nirvana.dat \ 
		-i HiSeq.10000.vcf \
		-o Hiseq_result_anno

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch nirvana.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. nirvana.swarm). For example:



```

nirvana -i vfc1.vcf -o vcf1 [... rest of the options]
nirvana -i vfc2.vcf -o vcf2 [... rest of the options]
nirvana -i vfc3.vcf -o vcf3 [... rest of the options]
nirvana -i vfc4.vcf -o vcf4 [... rest of the options]

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f nirvana.swarm --module nirvana
```

where


|  |  |
| --- | --- |
| --module nirvana Loads the nirvana module for each subjob in the swarm 
 | |








