

document.querySelector('title').textContent = 'pvactools on Biowulf';
pvactools on Biowulf


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



pVACtools is a cancer immunotherapy suite consisting of the following tools:

pVACseq
A cancer immunotherapy pipeline for identifying and prioritizing neoantigens from a list of tumor mutations.
pVACfuse
A tool for detecting neoantigens resulting from gene fusions.
pVACvector
A tool designed to aid specifically in the construction of DNA vector-based cancer vaccines.
pVACbind
A cancer immunotherapy pipeline for identifying and prioritizing neoantigens from a FASTA file.




### References:



Jasreet Hundal, Susanna Kiwala, Joshua McMichael, Christopher A Miller, Alexander T Wollam, Huiming Xia, Connor J Liu, Sidi Zhao, Yang-Yang Feng, Aaron P Graubert, Amber Z Wollam, Jonas Neichin, Megan Neveau, Jason Walker, William E Gillanders, Elaine R Mardis, Obi L Griffith, Malachi Griffith. [pVACtools: a computational toolkit to select and visualize cancer neoantigens](https://www.ncbi.nlm.nih.gov/pubmed/31907209). *Cancer Immunology Research*. 2020 Mar;8(3):409-420. DOI: 10.1158/2326-6066.


Documentation
* [pVACtools documentation](https://pvactools.readthedocs.io/en/latest/index.html) at readthedocs.io
* [pvactools website](http://www.pvactools.org)


Important Notes
* Module Name: pvactools (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded
* pvactools on Biowulf is installed inside a [Singularity container](singularity.html). For the user, this should not make any difference, as the pvactools commands will be run exactly as before. 
* The IEDB (Immune Epitope Database and Analysis Resource) MHC Class I and Class II archives are in /opt/iedb within the container. 
* For older versions of pvactools, the IEDB archives are also in /usr/local/apps/pvactools/IEDB on Biowulf.
* Note that the pVACviz utility cannot be used directly on Biowulf. Jobs must be submitted using sinteractive, sbatch, or swarm.
* We highly recommend [allocating lscratch](https://hpc.nih.gov/docs/userguide.html#local) for pvacseq. If possible, run your analysis in /lscratch/$SLURM\_JOB\_ID/ itself and move the output directory to /data when the run completes


The default module version as of March 2021 is 2.0 or greater. The newer versions [break compatibility](https://pvactools.readthedocs.io/en/latest/releases/2_0.html#breaking-changes) with previous versions. Please update your workflow or explicitly load older modules, for example module load pvactools/1.5.5



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --gres=lscratch:5**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load pvactools**
/lscratch/46116226 exists
[+] Loading singularity  3.7.1  on cn3144
[+] Loading pvactools  2.0.1

[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn3144 ~]$ **pvacseq download\_example\_data .**

[user@cn3144 ~]$ **pvacseq run \
 pvacseq\_example\_data/input.vcf \
 Test \
 HLA-A\*02:01,HLA-B\*35:01,DRB1\*11:01 \
 MHCflurry MHCnuggetsI MHCnuggetsII NNalign NetMHC PickPocket SMM SMMPMBEC SMMalign \
 pvacseq\_example\_output \
 -e1 8,9,10 \
 -e2 15**
Executing MHC Class I predictions
Converting .vcf to TSV
Completed
Converting VCF to TSV
Completed
Generating Variant Peptide FASTA and Key File
Completed
Parsing the Variant Peptide FASTA and Key File
Completed
Calculating Manufacturability Metrics
Completed
Splitting TSV into smaller chunks
Splitting TSV into smaller chunks - Entries 1-24
Completed
Generating Variant Peptide FASTA and Key Files
Generating Variant Peptide FASTA and Key Files - Epitope Length 8 - Entries 1-48
Generating Variant Peptide FASTA and Key Files - Epitope Length 9 - Entries 1-48
Generating Variant Peptide FASTA and Key Files - Epitope Length 10 - Entries 1-48
Completed
[...]
Creating combined reports
Creating aggregated report
Completed
Running Binding Filters
Completed
Running Coverage Filters
Completed
Running Transcript Support Level Filter
Complete
Running Top Score Filter
Completed

Done: Pipeline finished successfully. File /lscratch/46116226//pvacseq_example_output/combined/Test.filtered.tsv contains list of filtered putative neoantigens

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. pvactools.sh). For example:



```

#! /bin/bash
set -e

cd /data/$USER/somedir

module load pvactools
allele="HLA-A*23:01,HLA-A*68:02,HLA-B*07:17,HLA-B*08:01,HLA-C*02:02,HLA-C*17:01"
algorithms="MHCflurry MHCnuggetsI MHCnuggetsII NNalign NetMHC PickPocket SMM SMMPMBEC SMMalign"
pvacseq run input.vcf \
        Test \
        $allele \
        $algorithms \
        output_dir \
        -e1 8,9,10 \
        -e2 15

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#G] [--time=##:##:##] [--gres=lscratch:#] pvactools.sh
```


Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. pvactools.swarm). For example:



```

pvacseq run -e1 8,9,10 -e2 15 input1.vcf Test1 HLA-C*17:01 MHCflurry MHCnuggetsI MHCnuggetsII output1
pvacseq run -e1 8,9,10 -e2 15 input2.vcf Test2 HLA-C*17:01 MHCflurry MHCnuggetsI MHCnuggetsII output2
pvacseq run -e1 8,9,10 -e2 15 input3.vcf Test3 HLA-C*17:01 MHCflurry MHCnuggetsI MHCnuggetsII output3
pvacseq run -e1 8,9,10 -e2 15 input4.vcf Test4 HLA-C*17:01 MHCflurry MHCnuggetsI MHCnuggetsII output4

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f pvactools.swarm [-g #] --module pvactools --export=
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --module pvactools Loads the pvactools module for each subjob in the swarm 
 | |
 | |










