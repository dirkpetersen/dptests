

document.querySelector('title').textContent = 'stripenn: detection of atchitectural stripes from chromatin conformation capture (3C) data.';
**stripenn: detection of atchitectural stripes from chromatin conformation capture (3C) data.**


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



Stripenn is a command line interface python package developed for detection of atchitectural stripes 
from chromatin conformation capture (3C) data. It implements an algorithm rooted in computer vision 
for demarcation and quantification of the architectural stripes. Stripenn was demonstrated to outperform 
existing methods, be applicable in the context of analysis of B and T lymphocytes, and to allow 
examination of the role of sequence variation on the architectural stripes by studying the conservation 
of these features in inbred strains of mice. 



### References:


 Sora Yoon, Golnaz Vahedi   

*Stripenn detects architectural stripes from chromatin conformation data using computer vision*    

[bioRxiv (2021)](https://www.biorxiv.org/content/10.1101/2021.04.16.440239v1.full), doi: https://doi.org/10.1101/2021.04.16.440239.


Documentation
* [stripenn Github page](https://github.com/ysora/stripenn)


Important Notes
* Module Name: stripenn (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **STRIPENN\_HOME**  installation directory
	+ **STRIPENN\_BIN**       executable directory
	+ **STRIPENN\_DATA**    sample data folder



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. 
  
Sample session on a GPU node:



```

[user@biowulf ~]$ **sinteractive --mem=48g -c50 --gres=lscratch:20**
[user@cn2379 ~]$ **module load stripenn**
+] Loading singularity  3.8.4  on cn2379
[+] Loading stripenn  1.1.50

```

Basic usage:

```

[user@cn2379 ~]$ **stripenn --help**
Usage: stripenn [OPTIONS] COMMAND [ARGS]...

Options:
  --install-completion [bash|zsh|fish|powershell|pwsh]
                                  Install completion for the specified shell.
  --show-completion [bash|zsh|fish|powershell|pwsh]
                                  Show completion for the specified shell, to
                                  copy it or customize the installation.
  --help                          Show this message and exit.

Commands:
  compute   Finds stripe coordinates from 3D genomic data
  score     Calculates p-value and stripiness of given stripes based on...
  seeimage  Draws heatmap image of given position and color saturation...
[user@cn2379 ~]$ **stripenn compute --help**
Usage: stripenn compute [OPTIONS]

  Finds stripe coordinates from 3D genomic data

Options:
  --cool TEXT             Path to cool file  [required]
  -o, --out TEXT          Path to output directory  [required]
  --norm TEXT             Normalization method. It should be one of the column
                          name of Cooler.bin(). Check it with
                          Cooler.bins().columns (e.g., KR, VC, VC_SQRT)
                          [default: KR]
  -k, --chrom TEXT        Set of chromosomes. e.g., 'chr1,chr2,chr3', 'all'
                          will generate stripes from all chromosomes
                          [default: all]
  -c, --canny FLOAT       Canny edge detection parameter.  [default: 2.5]
  -l, --minL INTEGER      Minimum length of stripe.  [default: 10]
  -w, --maxW INTEGER      Maximum width of stripe.  [default: 8]
  -m, --maxpixel TEXT     Percentiles of the contact frequency data to
                          saturate the image. Separated by comma  [default:
                          0.95,0.96,0.97,0.98,0.99]
  -n, --numcores INTEGER  The number of cores will be used.  [default: 72]
  -p, --pvalue FLOAT      P-value cutoff for stripe.  [default: 0.1]
  --mask TEXT             Column coordinates to be masked. e.g.,
                          chr9:12345678-12345789  [default: 0]
  -s                      Use if system memory is low.
  --help                  Show this message and exit.

```

Download sample data:

```

[user@cn2379 ~]$ **wget https://www.dropbox.com/s/1bb2npvrzp3by5y/BL6.DPT.chr16.mcool?dl=0 -O test.mcool --no-check-certificate**

```

Run stripenn on the sample data:

```

[user@cn2379 ~]$ **stripenn compute --cool test.mcool::resolutions/25000 --out output\_dir/ -n 50 -k 16 -m 0.95,0.96,0.97,0.98,0.99**
Result will be stored in output_dir/
#######################################
Maximum pixel value calculation ...
#######################################
#######################################
Expected value calculation ...
#######################################

Processing Chromosome: 16
100%|███████████████████████████████████████████████████████████████████| 10/10 [00:00<00:00, 501.33it/s]
#######################################
Background distribution estimation ...
#######################################

1. Calculating the number of available columns ...

100%|████████████████████████████████████████████████████████████████████| 1/1 [00:00<00:00, 1111.07it/s]
2. Constituting background ...

100%|████████████████████████████████████████████████████████████████████| 1/1 [00:00<00:00, 1625.70it/s]
/opt/conda/envs/stripenn/lib/python3.9/site-packages/numpy/core/fromnumeric.py:3440: RuntimeWarning: Mean of empty slice.
  return _methods._mean(a, axis=axis, dtype=dtype,
Elapsed time for background estimation: 0.552 min
#################################################
Finding candidate stripes from each chromosome ...
#################################################

Chromosome: 16 / Maximum pixel: 95.0%
100%|█████████████████████████████████████████████████████████████████| 20/20 [00:00<00:00, 29610.34it/s]
Chromosome: 16 / Maximum pixel: 96.0%
100%|█████████████████████████████████████████████████████████████████| 20/20 [00:00<00:00, 15563.28it/s]
Chromosome: 16 / Maximum pixel: 97.0%
100%|█████████████████████████████████████████████████████████████████| 20/20 [00:00<00:00, 22795.13it/s]
Chromosome: 16 / Maximum pixel: 98.0%
100%|█████████████████████████████████████████████████████████████████| 20/20 [00:00<00:00, 22233.26it/s]
Chromosome: 16 / Maximum pixel: 99.0%
100%|█████████████████████████████████████████████████████████████████| 20/20 [00:00<00:00, 33367.57it/s]
##########################
Stripiness calculation ...
##########################

16
  0%|                                                                            | 0/244 [00:00<?, ?it/s]/opt/conda/envs/stripenn/lib/python3.9/site-packages/numpy/core/fromnumeric.py:3440: RuntimeWarning: Mean of empty slice.
  return _methods._mean(a, axis=axis, dtype=dtype,
/opt/conda/envs/stripenn/lib/python3.9/site-packages/numpy/core/_methods.py:189: RuntimeWarning: invalid value encountered in double_scalars
  ret = ret.dtype.type(ret / rcount)
100%|█████████████████████████████████████████████████████████████████| 244/244 [00:00<00:00, 313.91it/s]

1.249min taken.
Check the result stored in output_dir/
[user@cn2379 ~]$ **ls -l output\_dir**
-rw-r--r-- 1 user staff   180 Nov 23 10:52 result_filtered.txt
-rw-r--r-- 1 user staff 26578 Nov 23 10:52 result_unfiltered.txt
-rw-r--r-- 1 user staff   178 Nov 23 10:51 stripenn.log

```

End the interactive session:

```

[user@cn2379 ~]$ **exit**

```

Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. stripenn.sh). For example:



```

#!/bin/bash
set -e
module load stripenn
stripenn compute --cool test.mcool::resolutions/25000 --out output_dir/ -n 50 -k 16 -m 0.95,0.96,0.97,0.98,0.99

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
**sbatch stripenn.sh**
```





