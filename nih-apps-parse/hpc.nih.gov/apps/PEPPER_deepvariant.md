

document.querySelector('title').textContent = 'PEPPER\_deepvariant: long-read variant calling and nanopore assembly polishing';
**PEPPER\_deepvariant: long-read variant calling and nanopore assembly polishing**


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



PEPPER-Margin-DeepVariant is a haplotype-aware variant calling pipeline for processing third-generation 
nanopore sequence data. It outperforms the short-read-based single-nucleotide-variant identification method 
at the whole-genome scale and produces high-quality single-nucleotide variants in segmental duplications 
and low-mappability regions where short-read-based genotyping fails. It also can provide highly contiguous 
phase blocks across the genome with nanopore reads and perform de novo assembly polishing to produce 
diploid assemblies with high accuracy. This pipeline is also applicable to PacBio HiFi data, 
providing an efficient solution with superior performance over the current WhatsHap-DeepVariant standard. 



### References:


* Kishwar Shafin, Trevor Pesout, Pi-Chuan Chang, Maria Nattestad, Alexey Kolesnikov,
Sidharth Goel, Gunjan Baid, Mikhail Kolmogorov, Jordan M. Eizenga, Karen H. Miga,
Paolo Carnevali, Miten Jain, Andrew Carroll and Benedict Paten,   

*Haplotype-aware variant calling with PEPPER-Margin-DeepVariant enables high
accuracy in nanopore long-reads.*   

[Nature Methods volume 18, pages 1322–1332 (2021).](https://www.nature.com/articles/s41592-021-01299-w)


Documentation
* [PEPPER Github page](https://github.com/kishwarshafin/pepper)
* [PEPPER-Margin-DeepVariant r0.7 method description](https://github.com/kishwarshafin/pepper/blob/r0.7/docs/misc/pepper_v0.7_method_update.md)
* [PEPPER-Margin-DeepVariant usage and parameters](https://github.com/kishwarshafin/pepper/blob/r0.7/docs/usage/usage_and_parameters.md)


Important Notes
* Module Name: PEPPER\_deepvariant (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Implemented as a Singularity container
* Unusual environment variables set
	+ **PEPPER\_HOME**  installation directory
	+ **PEPPER\_BIN**    executables directory
	+ **PEPPER\_CONFIG**sample configuration files directory
	+ **PEPPER\_SRC**    source directory
	+ **PEPPER\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=8g -c10 --gres=lscratch:20**
[user@cn3316 ~]$ **module load PEPPER\_deepvariant** 
[+] Loading singularity  3.8.5-1  on cn4235
[+] Loading pepper_deepvariant 0.7  ...

```

The basic usage command for PEPPER\_deepvariant is:

```

[user@cn3123 user]$ **run\_pepper\_margin\_deepvariant -h**
usage: run_pepper_margin_deepvariant [-h] [--version] {call_variant} ...

Run PEPPER-Margin-DeepVariant for variant calling.Example run: run_pepper_margin_deepvariant -b  -f  -o  -t  <--ont\_r9\_guppy5\_sup/--ont\_r10\_q20/--hifi>

positional arguments:
 {call\_variant}

optional arguments:
 -h, --help show this help message and exit
 --version Show version.

[user@cn3123 user]$ **run\_pepper\_margin\_deepvariant call\_variant -h**
usage: run\_pepper\_margin\_deepvariant call\_variant -b BAM -f FASTA -o OUTPUT\_DIR -t THREADS [-r REGION] [-g] [--gvcf]
 [-s SAMPLE\_NAME] [--pepper\_snp\_skip\_indels]
 [--no\_pepper\_snp\_skip\_indels] [--use\_pepper\_hp] [--no\_use\_pepper\_hp]
 [--haplotag\_with\_pepper\_hp] [--no\_haplotag\_with\_pepper\_hp]
 [--use\_pepper\_snp\_candidates] [--no\_use\_pepper\_snp\_candidates]
 [--only\_pepper] [--only\_haplotag] [--phased\_output]
 [--skip\_final\_phased\_bam] [-p OUTPUT\_PREFIX] [-k] [--dry]
 [--pepper\_model PEPPER\_MODEL] [--pepper\_hp\_model PEPPER\_HP\_MODEL]
 [--pepper\_quantized] [--no\_pepper\_quantized]
 [--pepper\_downsample\_rate PEPPER\_DOWNSAMPLE\_RATE]
 [--pepper\_region\_size PEPPER\_REGION\_SIZE]
 [--pepper\_include\_supplementary] [--pepper\_min\_mapq PEPPER\_MIN\_MAPQ]
 [--pepper\_min\_snp\_baseq PEPPER\_MIN\_SNP\_BASEQ]
 [--pepper\_min\_indel\_baseq PEPPER\_MIN\_INDEL\_BASEQ]
 [--pepper\_snp\_frequency PEPPER\_SNP\_FREQUENCY]
 [--pepper\_insert\_frequency PEPPER\_INSERT\_FREQUENCY]
 [--pepper\_delete\_frequency PEPPER\_DELETE\_FREQUENCY]
 [--pepper\_min\_coverage\_threshold PEPPER\_MIN\_COVERAGE\_THRESHOLD]
 [--pepper\_candidate\_support\_threshold PEPPER\_CANDIDATE\_SUPPORT\_THRESHOLD]
 [--pepper\_snp\_candidate\_frequency\_threshold PEPPER\_SNP\_CANDIDATE\_FREQUENCY\_THRESHOLD]
 [--pepper\_indel\_candidate\_frequency\_threshold PEPPER\_INDEL\_CANDIDATE\_FREQUENCY\_THRESHOLD]
 [--pepper\_skip\_indels]
 [--pepper\_allowed\_multiallelics PEPPER\_ALLOWED\_MULTIALLELICS]
 [--pepper\_snp\_p\_value PEPPER\_SNP\_P\_VALUE]
 [--pepper\_insert\_p\_value PEPPER\_INSERT\_P\_VALUE]
 [--pepper\_delete\_p\_value PEPPER\_DELETE\_P\_VALUE]
 [--pepper\_snp\_q\_cutoff PEPPER\_SNP\_Q\_CUTOFF]
 [--pepper\_indel\_q\_cutoff PEPPER\_INDEL\_Q\_CUTOFF]
 [--pepper\_report\_snp\_above\_freq PEPPER\_REPORT\_SNP\_ABOVE\_FREQ]
 [--pepper\_report\_indel\_above\_freq PEPPER\_REPORT\_INDEL\_ABOVE\_FREQ]
 [--margin\_haplotag\_model MARGIN\_HAPLOTAG\_MODEL]
 [--margin\_hp\_haplotag\_model MARGIN\_HP\_HAPLOTAG\_MODEL]
 [--margin\_phase\_model MARGIN\_PHASE\_MODEL] [--dv\_model DV\_MODEL]
 [--dv\_alt\_aligned\_pileup DV\_ALT\_ALIGNED\_PILEUP]
 [--dv\_realign\_reads DV\_REALIGN\_READS]
 [--dv\_min\_mapping\_quality DV\_MIN\_MAPPING\_QUALITY]
 [--dv\_min\_base\_quality DV\_MIN\_BASE\_QUALITY]
 [--dv\_sort\_by\_haplotypes DV\_SORT\_BY\_HAPLOTYPES]
 [--dv\_parse\_sam\_aux\_fields DV\_PARSE\_SAM\_AUX\_FIELDS]
 [--dv\_add\_hp\_channel DV\_ADD\_HP\_CHANNEL]
 [--dv\_use\_hp\_information DV\_USE\_HP\_INFORMATION]
 [--dv\_use\_multiallelic\_mode DV\_USE\_MULTIALLELIC\_MODE]
 (--ont\_r9\_guppy5\_sup | --ont\_r10\_q20 | --hifi) [-h]

... 

[user@cn3123 user]$ **cp -r $PEPPER\_DATA/\* .**
[user@cn3123 user]$ **run\_pepper\_margin\_deepvariant call\_variant -b NA12878\_S1.chr20.10\_10p1mb.bam -f ucsc.hg19.chr20.unittest.fasta -o output -t 2 --hifi**
[08-16-2022 14:08:55] INFO: VARIANT CALLING MODULE SELECTED
[08-16-2022 14:08:55] INFO: [1/8] RUNNING THE FOLLOWING COMMAND
-------
mkdir -p output;
mkdir -p output/logs;
mkdir -p output/intermediate\_files;
cp /opt/pepper\_models/PEPPER\_VARIANT\_HIFI\_V7.pkl output/intermediate\_files/
-------
[08-16-2022 14:08:57] INFO: [2/8] RUNNING THE FOLLOWING COMMAND
-------
time pepper\_variant call\_variant -b NA12878\_S1.chr20.10\_10p1mb.bam -f ucsc.hg19.chr20.unittest.fasta -t 2 -m output/intermediate\_files/PEPPER\_VARIANT\_HIFI\_V7.pkl -o output/pepper/ --no\_quantized -s Sample --hifi 2>&1 | tee output/logs/1\_pepper.log
-------
[08-16-2022 14:08:59] INFO: HiFi VARIANT CALLING MODE SELECTED.
[08-16-2022 14:08:59] INFO: MODE: PEPPER SNP
[08-16-2022 14:08:59] INFO: THRESHOLDS ARE SET TO:
[08-16-2022 14:08:59] INFO: MIN MAPQ: 5
[08-16-2022 14:08:59] INFO: MIN SNP BASEQ: 10
[08-16-2022 14:08:59] INFO: MIN INDEL BASEQ: 10
[08-16-2022 14:08:59] INFO: MIN SNP FREQUENCY: 0.1
[08-16-2022 14:08:59] INFO: MIN INSERT FREQUENCY: 0.12
[08-16-2022 14:08:59] INFO: MIN DELETE FREQUENCY: 0.12
[08-16-2022 14:08:59] INFO: MIN COVERAGE THRESHOLD: 2
[08-16-2022 14:08:59] INFO: MIN CANDIDATE SUPPORT: 2
[08-16-2022 14:08:59] INFO: MIN SNP CANDIDATE FREQUENCY: 0.1
[08-16-2022 14:08:59] INFO: MIN INDEL CANDIDATE FREQUENCY: 0.12
[08-16-2022 14:08:59] INFO: SKIP INDEL CANDIDATES: False
[08-16-2022 14:08:59] INFO: MAX ALLOWED CANDIDATE IN ONE SITE: 4
[08-16-2022 14:08:59] INFO: MIN SNP PREDICTIVE VALUE: 0.01
[08-16-2022 14:08:59] INFO: MIN INSERT PREDICTIVE VALUE: 0.001
[08-16-2022 14:08:59] INFO: MIN DELETE PREDICTIVE VALUE: 0.01
[08-16-2022 14:08:59] INFO: SNP QV CUTOFF FOR RE-GENOTYPING: 20
[08-16-2022 14:08:59] INFO: INDEL QV CUTOFF FOR RE-GENOTYPING: 50
[08-16-2022 14:08:59] INFO: REPORT ALL SNPs ABOVE THRESHOLD: 0
[08-16-2022 14:08:59] INFO: REPORT ALL INDELs ABOVE THRESHOLD: 0
[08-16-2022 14:08:59] INFO: CALL VARIANT MODULE SELECTED
[08-16-2022 14:08:59] INFO: RUN-ID: 08162022\_140859
[08-16-2022 14:08:59] INFO: IMAGE OUTPUT: output/pepper/images\_08162022\_140859/
[08-16-2022 14:08:59] INFO: STEP 1/3 GENERATING IMAGES:
[08-16-2022 14:08:59] INFO: COMMON CONTIGS FOUND: ['chr20']
[08-16-2022 14:08:59] INFO: TOTAL CONTIGS: 1 TOTAL INTERVALS: 631 TOTAL BASES: 63025519
[08-16-2022 14:08:59] INFO: STARTING PROCESS: 0 FOR 316 INTERVALS
[08-16-2022 14:08:59] INFO: [THREAD 00] 10/316 COMPLETE (3%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 20/316 COMPLETE (6%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 30/316 COMPLETE (9%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 40/316 COMPLETE (12%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 50/316 COMPLETE (15%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 60/316 COMPLETE (18%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 70/316 COMPLETE (22%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 80/316 COMPLETE (25%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 90/316 COMPLETE (28%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 100/316 COMPLETE (31%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 110/316 COMPLETE (34%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 120/316 COMPLETE (37%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 130/316 COMPLETE (41%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 140/316 COMPLETE (44%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 150/316 COMPLETE (47%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 160/316 COMPLETE (50%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 170/316 COMPLETE (53%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 180/316 COMPLETE (56%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 190/316 COMPLETE (60%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 200/316 COMPLETE (63%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 210/316 COMPLETE (66%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 220/316 COMPLETE (69%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 230/316 COMPLETE (72%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 240/316 COMPLETE (75%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 250/316 COMPLETE (79%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 260/316 COMPLETE (82%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 270/316 COMPLETE (85%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:08:59] INFO: [THREAD 00] 280/316 COMPLETE (88%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:09:00] INFO: [THREAD 00] 290/316 COMPLETE (91%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:09:00] INFO: [THREAD 00] 300/316 COMPLETE (94%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:09:00] INFO: [THREAD 00] 310/316 COMPLETE (98%) [ELAPSED TIME: 0 Min 0 Sec]
[08-16-2022 14:09:00] INFO: THREAD 0 FINISHED SUCCESSFULLY.
[08-16-2022 14:09:00] INFO: FINISHED IMAGE GENERATION
[08-16-2022 14:09:00] INFO: TOTAL ELAPSED TIME FOR GENERATING IMAGES: 0 Min 0 Sec
[08-16-2022 14:09:00] INFO: STEP 2/3 RUNNING INFERENCE
[08-16-2022 14:09:00] INFO: OUTPUT: output/pepper/predictions\_08162022\_140859/
[08-16-2022 14:09:00] INFO: DISTRIBUTED CPU SETUP.
[08-16-2022 14:09:00] INFO: TOTAL CALLERS: 2
[08-16-2022 14:09:00] INFO: THREADS PER CALLER: 1
[08-16-2022 14:09:00] INFO: MODEL LOADING TO ONNX
[08-16-2022 14:09:00] INFO: SETTING THREADS TO: 1.
[08-16-2022 14:09:00] INFO: STARTING INFERENCE.
[08-16-2022 14:09:00] INFO: TOTAL SUMMARIES: 2.
[08-16-2022 14:09:00] INFO: SUMMARY PROCESSED 1/2.
[08-16-2022 14:09:00] INFO: SUMMARY PROCESSED 2/2.
[08-16-2022 14:09:00] INFO: THREAD 0 FINISHED SUCCESSFULLY.
[08-16-2022 14:09:00] INFO: FINISHED PREDICTION
[08-16-2022 14:09:00] INFO: ELAPSED TIME: 0 Min 0 Sec
[08-16-2022 14:09:00] INFO: PREDICTION FINISHED SUCCESSFULLY.
[08-16-2022 14:09:00] INFO: TOTAL ELAPSED TIME FOR INFERENCE: 0 Min 0 Sec
[08-16-2022 14:09:00] INFO: STEP 3/3 FINDING CANDIDATES
[08-16-2022 14:09:00] INFO: OUTPUT: output/pepper/
[08-16-2022 14:09:00] INFO: STARTING CANDIDATE FINDING.
[08-16-2022 14:09:01] INFO: FINISHED PROCESSING, TOTAL CANDIDATES FOUND: 194
[08-16-2022 14:09:01] INFO: FINISHED PROCESSING, TOTAL VARIANTS IN PEPPER: 16
[08-16-2022 14:09:01] INFO: FINISHED PROCESSING, TOTAL VARIANTS SELECTED FOR RE-GENOTYPING: 178
[08-16-2022 14:09:01] INFO: TOTAL TIME SPENT ON CANDIDATE FINDING: 0 Min 0 Sec
[08-16-2022 14:09:01] INFO: TOTAL ELAPSED TIME FOR FINDING CANDIDATES: 0 Min 1 Sec

real 0m4.129s
user 0m2.901s
sys 0m1.906s
[08-16-2022 14:09:01] INFO: [3/8] RUNNING THE FOLLOWING COMMAND
-------
mv output/pepper/PEPPER\_VARIANT\_FULL.vcf output/;
mv output/pepper/PEPPER\_VARIANT\_OUTPUT\_PEPPER.vcf output/;
mv output/pepper/PEPPER\_VARIANT\_OUTPUT\_VARIANT\_CALLING.vcf output/;
bgzip output/PEPPER\_VARIANT\_FULL.vcf;
tabix -p vcf output/PEPPER\_VARIANT\_FULL.vcf.gz;
bgzip output/PEPPER\_VARIANT\_OUTPUT\_PEPPER.vcf;
tabix -p vcf output/PEPPER\_VARIANT\_OUTPUT\_PEPPER.vcf.gz;
bgzip output/PEPPER\_VARIANT\_OUTPUT\_VARIANT\_CALLING.vcf;
tabix -p vcf output/PEPPER\_VARIANT\_OUTPUT\_VARIANT\_CALLING.vcf.gz;
rm -rf output/pepper/;
echo "CONTIGS FOUND IN PEPPER VCF:";
zcat output/PEPPER\_VARIANT\_FULL.vcf.gz | grep -v '#' | cut -f1 | uniq
...

```


```

[user@cn3316 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. fithichip.sh). For example:



```

#!/bin/bash
module load PEPPER_deepvariant

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] fithichip.sh
```





