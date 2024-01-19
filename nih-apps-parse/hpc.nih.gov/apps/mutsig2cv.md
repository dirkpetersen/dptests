

document.querySelector('title').textContent = 'MutSig2CV on Biowulf';
MutSig2CV on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



mutsig2cv analyzes somatic point mutations discovered in DNA sequencing, identifying genes mutated more often than expected by chance.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --mem=10g --gres=lscratch:100**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **cd /lscratch/${SLURM\_JOB\_ID}**
[user@cn3144 ~]$ **module load mutsig2cv**
[user@cn3144 ~]$ **cp /usr/local/apps/mutsig2cv/TEST\_DATA/MutSigCV\_example\_data.1.0.1.zip .**
[user@cn3144 ~]$ **unzip MutSigCV\_example\_data.1.0.1.zip**
[user@cn3144 ~]$ **cd MutSigCV\_example\_data.1.0.1**
[user@cn3144 ~]$ **cp /usr/local/matlab-compiler/v81.tar.gz .**
[user@cn3144 ~]$ **tar -xzf v81.tar.gz**
[user@cn3144 ~]$ **cp /usr/local/apps/mutsig2cv/3.11/reference.tar.gz .**
[user@cn3144 ~]$ **tar -xzf reference.tar.gz**
[user@cn3144 ~]$ **run\_MutSig2CV.sh v81 LUSC.mutations.maf TEST\_OUTPUT**
[...]
16:16:47  16:41:52  8187/18862 ITGA9       3       1000  CV 0.617212 CL 1.000000 FN 0.326000  [625  81   10   1    0    1]
16:16:47  16:41:53  8188/18862 ITGAD       9       1000  CV 0.470805 CL 1.000000 FN 0.349000  [625  81   10   1    0    1]
16:16:48  16:41:53  8189/18862 ITGAE       4       1000  CV 0.918407 CL 1.000000 FN 0.731000  [625  81   10   1    0    1]
16:16:48  16:41:53  8190/18862 ITGAL       1          0  CV 0.998619 CL 1.000000 FN 1.000000  [625  81   10   1    0    1]
16:16:48  16:41:53  8191/18862 ITGAM       7       1000  CV 0.653397 CL 1.000000 FN 0.068000  [625  81   10   1    0    1]
16:16:48  16:41:54  8192/18862 ITGAV       4       1000  CV 0.849691 CL 1.000000 FN 0.100000  [626  81   10   1    0    1]
16:16:48  16:41:54  8193/18862 ITGAX       8       1000  CV 0.410713 CL 1.000000 FN 0.885000  [627  81   10   1    0    1]
16:16:52  16:41:56  8210/18862 ITIH4       2       1000  CV 0.872721 CL 1.000000 FN 0.698000  [629  82   10   1    0    1]
16:16:52  16:41:56  8211/18862 ITIH5       5       1000  CV 0.717304 CL 1.000000 FN 0.729000  [629  82   10   1    0    1]
[...]
16:44:08  16:44:08 18860/18862 ZYX         1          0  CV 0.817634 CL 1.000000 FN 1.000000  [1591 184  19   4    1    7]
16:44:08  16:44:08 18861/18862 ZZEF1       2       1000  CV 0.980706 CL 1.000000 FN 0.500000  [1591 184  19   4    1    7]
16:44:08  16:44:08 18862/18862 ZZZ3        6       1000  CV 0.302995 CL 0.055000 FN 0.507000  [1591 184  19   4    1    7]

9 genes with q<=0.1
Saving results...   [save_struct] 38/70 39/70 40/70 41/70 42/70 43/70 44/70 45/70 46/70 47/70 48/70 49/70 50/70 51/70 52/70 
53/70 54/70 55/70 56/70 57/70 58/70 59/70 60/70 61/70 62/70 63/70 64/70 65/70 66/70 67/70 68/70 69/70 70/70  [collapse] [write]
Done.

[user@cn3144 ~]$ 
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. mutsig2cv.sh). For example:



```

#!/bin/bash
set -e
module load mutsig2cv
cd /lscratch/${SLURM_JOB_ID}
cp /usr/local/apps/mutsig2cv/TEST_DATA/MutSigCV_example_data.1.0.1.zip .
unzip MutSigCV_example_data.1.0.1.zip
cd MutSigCV_example_data.1.0.1
cp /usr/local/matlab-compiler/v81.tar.gz .
tar -xzf v81.tar.gz
cp /usr/local/apps/mutsig2cv/3.11/reference.tar.gz .
tar -xzf reference.tar.gz
run_MutSig2CV.sh v81 LUSC.mutations.maf TEST_OUTPUT

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --mem=10g --gres=lscratch:100 mutsig2cv.sh
```


10 GB memory is sufficient for this example job. You may need to increase the memory allocation for your own jobs. 



