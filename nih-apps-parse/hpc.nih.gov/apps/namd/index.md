

document.querySelector('title').textContent = ' NAMD 2.12 Benchmarks';

NAMD Benchmarks



|  |
| --- |
| 
Quick Links
[NAMD 3.0 on GPUs](#3.0-GPU)
[NAMD 2.14 on CPUs](#2.14-CPU)
[NAMD 2.14 on GPUs](#2.14-GPU)
 |


**Important:** Please read the webpage [Making efficient use of Biowulf's Multinode Partition](https://hpc.nih.gov/policies/multinode.html) before running large parallel jobs. 

Older Benchmarks: [[NAMD 2.12](2.12.html)] [[NAMD 2.10](2.10.html)]


NAMD 3.0 beta on Biowulf GPUS

**April 2023** STMV benchmark, 1,066,628 atoms, periodic, PME (available from [here](http://www.ks.uiuc.edu/Research/namd/utilities/stmv.tar.gz))

[NAMD config file used for these GPU benchmarks](new-stmv.namd.5ksteps) was 
run with these parameters: P100 & V100: GPUPERF=1; V100X & A100: GPUPERF=2. (as recommended by 
[David Hardy, UIUC](https://www.ks.uiuc.edu/~dhardy/)). 

These benchmarks are simply provided here to give an estimate of the performance. Users who use NAMD 3.0 for their own runs
should carefully read the [NAMD 3.0 
documentation](https://www.ks.uiuc.edu/Research/namd/3.0/announce.html).

Hardware:
* p100: Intel(R)\_Xeon(R)\_CPU\_E5-2680\_v4\_@\_2.40GHz, Tesla\_P100-PCIE-16GB 
* v100: Intel(R)\_Xeon(R)\_CPU\_E5-2680\_v4\_@\_2.40GHz, Tesla\_V100-PCIE-16GB 
* v100x: Intel(R)\_Xeon(R)\_Gold\_6140@\_2.30GHz, Tesla\_V100-SXM2-32GB 
* a100: AMD EPYC 7543P@\_2.40GHz, Tesla\_A100-SXM4-80GB


OS: Rocky Linux release 8.7

![](3.0beta3-stmv.png)


NAMD 2.14 on CPUs

The Infiniband CPU benchmarks were performed with NAMD 2.14, Linux-x86\_64-ibverbs , downloaded from [here](http://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=NAMD).  

Hardware: 28 x 2.3 GHz (Intel E5-2695v3) (Biowulf multinode partition, x2695 nodes)

![](2.14-x2695.png)

NAMD 2.14 on Biowulf GPUs

GPU benchmarks were performed with NAMD 2.14, Linux-x86\_64-ibverbs-smp-CUDA downloaded from [here](http://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=NAMD).

 Hardware: as above for the 3.0alpha13 benchmarks.  

 OS: Rocky Linux release 8.7

![](2.14-stmv.png)






























