

document.querySelector('title').textContent = ' SNP2HLA on Biowulf';
 SNP2HLA on Biowulf



|  |
| --- |
| 
Quick Links
[Interactive Job on Biowulf](#int)
[Single Batch Job on Biowulf](#serial)
[Swarm of Jobs](#swarm)
[Documentation](#doc)
 |


SNP2HLA is a tool to impute amino acid polymorphisms and single nucleotide polymorphisms in human luekocyte antigenes (HLA) within the major histocompatibility complex (MHC) region in chromosome 6.
 The unique feature of SNP2HLA is that it imputes not only the classical HLA alleles but also the amino acid sequences of those classical alleles, so that individual amino acid sites can be directly tested for association. This allows for facile amino-acid focused downstream analysis.
 SNP2HLA also provides a companion package, MakeReference. This is software that builds the reference panel that can be used for SNP2HLA. This package is used for the situation that the provided reference panel is inappropriate (e.g. different populations), and the user wants to build the reference panel by him(her)self, e.g. typing the HLA alleles in a subset of individuals.
 SNP2HLA is developed by Sherman Jia and Buhm Han in the labs of [Soumya Raychaudhuri](http://www.broadinstitute.org/personal/soumya/)and [Paul de Bakker](http://www.broadinstitute.org/personal/debakker/)at the [Brigham and Women's Hospital](http://www.broadinstitute.org/mpg/snp2hla/www.brighamandwomens.org/) and [Harvard Medical School](http://hms.harvard.edu/), and the Broad Institute. 


* Example files for SNP2HLA can be copied from*/usr/local/apps/SNP2HLA/1.0.3/SNP2HLA/Example*
* Example files for MakeReference can be copied from*/usr/local/apps/SNP2HLA/1.0.3/MakeReference/Example*


 Interactive job 

[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load snp2hla**

[user@cn3144 ~]$ **cd /data/$USER/**

[user@cn3144 ~]$ **mkdir snp2hla; cp -r /usr/local/apps/snp2HLA/1.0.3/SNP2HLA/Example ./snp2hla**

[user@cn3144 ~]$ **cd /data/$USER/snp2hla/Example**

[user@cn3144 ~]$ **SNP2HLA.csh 1958BC HM\_CEU\_REF 1958BC\_IMPUTED plink 2000 1000**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

**Please Note, the above example uses '2000' in the command which will be passed into the SNP2HLA.csh, and therefore -Xmx2000m. If user needs more memory for analysis, please modify the commands in the scripts accordingly.** 

Running a single batch job on Biowulf



 1. Create a script file similar to the 
 lines below.



```

#!/bin/bash

module load snp2hla
cd /data/$USER/Examples
SNP2HLA.csh 1958BC HM_CEU_REF 1958BC_IMPUTED plink 4000 1000
```

2. Submit the script on biowulf:



```
$ sbatch jobscript
```


 For more memory requirement (default 4gb), use --mem flag: 
 
```
$ sbatch --mem=10g jobscript
```

Please Note, the above example uses '4000' in the command which will be passed into the SNP2HLA.csh, and therefore -Xmx4000m. If user needs more memory for analysis, please modify the commands in the scripts accordingly.
 
Running a swarm of jobs on Biowulf
Setup a swarm command file:
 
```

  cd /data/$USER/dir1; SNP2HLA.csh 1958BC HM_CEU_REF 1958BC_IMPUTED plink 1500 1000
  cd /data/$USER/dir2; SNP2HLA.csh 1958BC HM_CEU_REF 1958BC_IMPUTED plink 1500 1000
  cd /data/$USER/dir3; SNP2HLA.csh 1958BC HM_CEU_REF 1958BC_IMPUTED plink 1500 1000
	[......]
```

Submit the swarm file:
 
```

  $ swarm -f swarmfile --module snp2hla
```

 
 -f: specify the swarmfile name   

 --module: set environmental variables for each command line in the file
To allocate more memory, use -g flag: 

```
  $ swarm -f swarmfile -g 4 --module snp2hla
```

-g: allocate more memory


 For more information regarding running swarm, see [**swarm.html**](http://hpc.cit.nih.gov/apps/swarm.html)
**Please Note, the above example uses '1500' in the command which will be passed into the SNP2HLA.csh, and therefore -Xmx1500m. If user needs more memory for analysis, please modify the commands in the scripts accordingly.**
Please Note, the above example uses '4000' in the command which will be passed into the SNP2HLA.csh, and therefore -Xmx4000m. If user needs more memory for analysis, please modify the commands in the scripts accordingly.
Documentation
<http://www.broadinstitute.org/mpg/snp2hla/>




































