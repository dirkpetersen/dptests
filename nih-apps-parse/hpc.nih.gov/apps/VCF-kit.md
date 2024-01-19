

document.querySelector('title').textContent = 'VCF-kit: assorted utilities for the variant call format ';
**VCF-kit: assorted utilities for the variant call format** 


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



VCF-kit is a collection of utility tools for processing and analyzing the VCF 
(variant call format) files, including primer generation for variant 
validation, dendrogram production,genotype imputation from sequence data 
in linkage studies, and additional tools to be used by statistical and 
population geneticists. 



### References:


* Daniel E. Cook and Erik C. Andersen   

*VCF-kit: assorted utilities for the variant
call format*   

[Bioinformatics](https://academic.oup.com/bioinformatics/article/33/10/1581/2908861)  2017, **33**(10), 1581–1582. doi: 10.1093/bioinformatics/btx011


Documentation
* [VCF-kit home page](https://vcf-kit.readthedocs.io/en/latest/)


Important Notes
* Module Name: VCF-kit (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **VK\_HOME**  VCF-kit installation directory
	+ **VK\_BIN**       VCF-kit executable directory
	+ **VK\_DOC**       VCF-kit documentation directory
	+ **VK\_DATA**    VCF-kit test data directory
	+ **VK\_SRC**    VCF-kit source code directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=4g** 
[user@cn3316 ~]$**module load VCF-kit**
+] Loading singularity  3.8.5-1  on cn4213
[+] Loading VCF-kit 0.2.9  ...
[user@cn3316 ~]$**vk -h**
usage:
  vk  [...]
 vk setup
 vk -h | --help
 vk --version

commands:
 calc
 call
 filter
 geno
 genome
 hmm
 phylo
 primer
 rename
 tajima
 vcf2tsv
[user@cn3316 ~]$**vk calc -h**
usage:
 vk calc sample\_hom\_gt 
 vk calc genotypes [--frequency] 
 vk calc spectrum 

Example

options:
 -h --help Show this screen.
 --version Show version.
[user@cn3316 ~]$**cp $VK\_DATA/\* .** 
[user@cn3316 ~]$**vk calc genotypes test.vcf.gz** 
n ref het alt mis
937 14 0 0 0
327 13 0 0 1
243 13 0 1 0
168 12 0 0 2
101 11 0 0 3
93 12 0 1 1
90 12 0 2 0
73 10 0 0 4
63 11 0 3 0
47 10 0 4 0
37 11 0 1 2
36 9 0 0 5
36 8 0 0 6
32 11 0 2 1
31 10 0 1 3
28 10 0 3 1
28 9 0 4 1
25 10 0 2 2
23 13 1 0 0
21 9 0 5 0
20 7 0 0 7
19 9 0 1 4
19 8 0 1 5
17 8 0 6 0
15 11 1 0 2
15 6 0 0 8
13 9 0 2 3
13 9 0 3 2
...
1 7 5 2 0
1 3 3 1 7
1 8 4 2 0
1 0 0 12 2
[user@cn3316 ~]$**vk calc genotypes QX1211.indels.vcf.gz** 
n ref het alt mis
74493 0 0 1 0
8567 0 1 0 0

```

End the interactive session:

```

[user@cnR3316 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





