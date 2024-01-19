

document.querySelector('title').textContent = 'cryoDRGN: Deep Reconstructing Generative Networks for cryo-EM heterogeneous reconstruction';
**cryoDRGN: Deep Reconstructing Generative Networks for cryo-EM heterogeneous reconstruction**


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



CryoDRGN is an algorithm that leverages the representation power of deep neural networks 
to directly reconstruct continuous distributions of 3D density maps 
and map per-particle heterogeneity of single-particle cryo-EM datasets. 
It contains interactive tools to visualize a dataset’s distribution of per-particle 
variability, generate density maps for exploratory analysis, extract particle
subsets for use with other tools and generate trajectories to visualize molecular motions. 

CryoDRGN is open-source software
freely available at http://cryodrgn.csail.mit.edu.



### References:


* Ellen D. Zhong, Tristan Bepler, Bonnie Berger, and Joseph H. Davis,   

 *CryoDRGN: reconstruction of heterogeneous cryo-EM structures using neural networks.*   

[Nature Methods,](https://www.nature.com/articles/s41592-020-01049-4) VOL 18 | FebruarY 2021 | 176–185   
* Ellen D. Zhong, Tristan Bepler, Joseph H. Davis, and Bonnie Berger,   

*Reconstructing continuous distributions of 3D protein structure from cryo-EM images.*   

 [ICLR 2020, Spotlight presentation.](https://arxiv.org/abs/1909.05215)


Documentation
* [CryoDRGN User Guide](https://ez-lab.gitbook.io/cryodrgn/)
* [CryoDRGN Tutorial](https://www.notion.so/cryoDRGN-tutorial-b932c021cb2c415282f182048bac16ff)
* [CryoDRGN Main site](http://cryodrgn.csail.mit.edu)
* [CryoDRGN on github](https://github.com/zhonge/cryodrgn)


Important Notes
* Module Name: cryoDRGN (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded
* Unusual environment variables set
	+ **CRYODRGN\_HOME**   cryoDRGN installation directory
	+ **CRYODRGN\_BIN**         cryoDRGN executable directory
	+ **CRYODRGN\_TEST**     cryoDRGN test data directory* Example files in **$CRYODRGN\_TEST**



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --gres=gpu:k80:1**
[user@cn4199 ~]$ **module load cryoDRGN**
[+] Loading singularity  3.8.5-1  on cn4199
[+] Loading cryoDRGN 1.1.0  ...

```

Copy test data to your current folder:

```

[user@user@cn4199 ~]$ **cp -r $CRYODRGN\_TEST/\* .**

```

Downsample test data:

```

[user@cn4199 ~]$ **cryodrgn -h**
usage: cryodrgn [-h] [--version]
                {preprocess,downsample,parse_pose_csparc,parse_pose_star,parse_ctf_csparc,parse_ctf_star,train_nn,backproject_voxel,train_vae,eval_vol,eval_images,analyze,analyze_landscape,analyze_landscape_full,pc_traversal,graph_traversal,view_config}
                ...

CryoDRGN neural network reconstruction

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

Choose a command:
  {preprocess,downsample,parse_pose_csparc,parse_pose_star,parse_ctf_csparc,parse_ctf_star,train_nn,backproject_voxel,train_vae,eval_vol,eval_images,analyze,analyze_landscape,analyze_landscape_full,pc_traversal,graph_traversal,view_config}
[user@cn4199 ~]$ **cryodrgn downsample -h**
usage: cryodrgn downsample [-h] -D D -o MRCS [-b B] [--is-vol] [--chunk CHUNK]
                           [--datadir DATADIR] [--max-threads MAX_THREADS]
                           mrcs

Downsample an image stack or volume by clipping fourier frequencies

positional arguments:
  mrcs                  Input particles or volume (.mrc, .mrcs, .star, or
                        .txt)

optional arguments:
  -h, --help            show this help message and exit
  -D D                  New box size in pixels, must be even
  -o MRCS               Output projection stack (.mrcs)
  -b B                  Batch size for processing images (default: 5000)
  --is-vol              Flag if input .mrc is a volume
  --chunk CHUNK         Chunksize (in # of images) to split particle stack
                        when saving
  --datadir DATADIR     Optionally provide path to input .mrcs if loading from
                        a .star or .cs file
  --max-threads MAX_THREADS
                        Maximum number of CPU cores for parallelization
                        (default: 16)
[user@cn4199 ~]$ **cryodrgn downsample data/toy\_projections.mrcs -D 24 -o particles.24.mrcs**
2022-09-08 12:24:55     Processing batch 0
2022-09-08 12:24:55     (1000, 24, 24)
2022-09-08 12:24:55     Saving /gpfs/gsfs7/users/user/cryoDRGN/particles.24.mrcs

```

Parse image poses from a consensus homogeneous reconstruction:

```

[user@cn4199 ~]$ **cryodrgn parse\_pose\_star -h**
usage: cryodrgn parse_pose_star [-h] -o PKL [-D D] [--Apix APIX] input

Parse image poses from RELION .star file

positional arguments:
  input        RELION .star file

optional arguments:
  -h, --help   show this help message and exit
  -o PKL       Output pose.pkl

Optionally provide missing image parameters:
  -D D         Box size of reconstruction (pixels)
  --Apix APIX  Pixel size (A); Required if translations are specified in
               Angstroms
[user@cn4199 ~]$ **cryodrgn parse\_pose\_star data/relion31.star -o test.pkl -D 300 --Apix 1.03**
2022-09-08 12:26:14     5 particles
2022-09-08 12:26:14     Euler angles (Rot, Tilt, Psi):
2022-09-08 12:26:14     [-102.30296    82.318041  124.706463]
2022-09-08 12:26:14     Converting to rotation matrix:
2022-09-08 12:26:14     [[ 0.81941806 -0.10080704  0.56426234]
 [-0.53288075  0.22868946  0.81470193]
 [-0.21116853 -0.96826601  0.13367414]]
2022-09-08 12:26:14     Translations (pixels):
2022-09-08 12:26:14     [ 0.00688667 -0.23052744]
2022-09-08 12:26:14     Writing /gpfs/gsfs7/users/user/cryoDRGN/test.pkl

```

Parse CTF parameters from a .star/.cs file:

```

[user@cn4199 ~]$ **cryodrgn parse\_ctf\_star -h**
usage: cryodrgn parse_ctf_star [-h] -o O [--png PNG] [-D D] [--Apix APIX]
                               [--kv KV] [--cs CS] [-w W] [--ps PS]
                               star

Parse CTF parameters from a RELION .star file

positional arguments:
  star         Input

optional arguments:
  -h, --help   show this help message and exit
  -o O         Output pkl of CTF parameters
  --png PNG    Optionally plot the CTF

Optionally provide missing image parameters:
  -D D         Image size in pixels
  --Apix APIX  Angstroms per pixel
  --kv KV      Accelerating voltage (kV)
  --cs CS      Spherical abberation (mm)
  -w W         Amplitude contrast ratio
  --ps PS      Phase shift (deg)
[user@cn4199 ~]$ **cryodrgn parse\_ctf\_star data/relion31.star -D 300 --Apix 1.03 -o ctf.pkl --relio
n31 --kv 10 --cs 1 -w 0.5**
2021-02-24 10:24:39     5 particles
2021-02-24 10:24:39     Overriding accerlating voltage with 10.0 kV
2021-02-24 10:24:39     Overriding spherical abberation with 1.0 mm
2021-02-24 10:24:39     Overriding amplitude contrast ratio with 0.5
2021-02-24 10:24:39     CTF parameters for first particle:
2021-02-24 10:24:39     Image size (pix)  : 300
2021-02-24 10:24:39     A/pix             : 1.03
2021-02-24 10:24:39     DefocusU (A)      : 13108.082418
2021-02-24 10:24:39     DefocusV (A)      : 12845.582418
2021-02-24 10:24:39     Dfang (deg)       : -160.39
2021-02-24 10:24:39     voltage (kV)      : 10.0
2021-02-24 10:24:39     cs (mm)           : 1.0
2021-02-24 10:24:39     w                 : 0.5
2021-02-24 10:24:39     Phase shift (deg) : 0.0
2021-02-24 10:24:39     Saving ctf.pkl
[user@cn4199 ~]$  **cryodrgn downsample data/relion31.star -D 24 -o relion31.24.mrcs**
2022-09-08 12:33:30     Processing batch 0
2022-09-08 12:33:30     (5, 24, 24)
2022-09-08 12:33:30     Saving /gpfs/gsfs7/users/user/cryoDRGN/relion31.24.mrcs
[user@cn4199 ~]$ **exit**
salloc.exe: Relinquishing job allocation 59748321
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. CryoDRGN.sh). For example:



```

#!/bin/bash
module load CryoDRGN      
cp -r $CRYODRGN_TEST/* .
cryodrgn parse_pose_star data/relion31.star -o test.pkl -D 300  --Apix 1.03
cryodrgn parse_ctf_star data/relion31.star -D 300 --Apix 1.03 -o ctf.pkl --kv 10 --cs 1 -w 0.5

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] CryoDRGN.sh
```







