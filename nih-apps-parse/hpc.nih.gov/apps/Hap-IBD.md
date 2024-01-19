

document.querySelector('title').textContent = 'Hap-IBD: detecting identity-by-descent (IBD) segments and homozygosity-by-descent (HBD) segments';
**Hap-IBD: detecting identity-by-descent (IBD) segments and homozygosity-by-descent (HBD) segments**


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



The hap-ibd program detects identity-by-descent (IBD) segments and homozygosity-by-descent (HBD) segments 
in phased genotype data. The hap-ibd program can analyze data sets with hundreds of thousands of samples.



Documentation
* [hap-ibd GitHub page](https://github.com/browning-lab/hap-ibd)


Important Notes
* Module Name: hap-ibd (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **HAPIBD\_BIN**       executable directory
	+ **HAPIBD\_DATA**  sample data directory
	 + **JARPATH**  path to the folder containing JAR file



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
[user@cn3335 ~]$ **module load Hap-IBD**
[+] Loading java 12.0.1  ...
[+] Loading Hap-IBD  20221201
[user@cn3335 ~]$ **hap-ibd**
hap-ibd.jar  [ version 1.0, 20May22.818 ]

Syntax: java -jar hap-ibd.jar [arguments in format: parameter=value]

Data Parameters:
  gt=<VCF file with GT field>                         (required)
  map=<PLINK map file with cM units>                  (required)
  out=<output file prefix>                            (required)
  excludesamples=<excluded samples file>              (optional)

Algorithm Parameters:
  min-seed=<min cM length of seed segment>            (default: 2.0)
  max-gap=<max base pairs in non-IBS gap>             (default: 1000)
  min-extend=<min cM length of extension segment>     (default: min(1.0, min-seed))
  min-output=<min cM length of output segment>        (default: 2.0)
  min-markers=<min markers in seed segment>           (default: 100)
  min-mac=<minimum minor allele count filter>         (default: 2)
  nthreads=<number of computational threads>          (default: all CPU cores)
[user@cn3335 ~]$ **cp $HAPIBD\_DATA/\* .**
[user@cn3335 ~]$ **hap-ibd gt=target.truth.vcf.gz map=target.map out=hap-ibd.out**
Copyright (C) 2019 Brian L. Browning
Enter "java -jar hap-ibd.jar" to print a list of command line arguments

Program            :  hap-ibd.jar  [ version 1.0, 20May22.818 ]
Start Time         :  11:21 AM EST on 31 Jan 2023
Max Memory         :  30688 MB

Parameters
  gt               :  target.truth.vcf.gz
  map              :  target.map
  out              :  hap-ibd.out
  min-seed         :  2.0
  max-gap          :  1000
  min-extend       :  1.0
  min-output       :  2.0
  min-markers      :  100
  min-mac          :  2
  nthreads         :  24

Statistics
  samples          :  300
  markers          :  24965
  IBD segments     :  181
  IBD segs/sample  :  0.6
  HBD segments     :  0
  HBD segs/sample  :  0.000

Wallclock Time:    :  2 seconds
End Time           :  11:21 AM EST on 31 Jan 2023

```





