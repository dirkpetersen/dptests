

document.querySelector('title').textContent = 'CHESS: comparison and automatic feature extraction for chromatin contact data';
**CHESS: comparison and automatic feature extraction for chromatin contact data**


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



The CHESS (Comparison of Hi-C Experiments using Structural Similarity) application 
implements an algorithm for the comparison of chromatin contact maps and automatic
differential feature extraction. 



### References:


Silvia Galan, Nick Machnik, Kai Kruse, Noelia Díaz, Marc A. Marti-Renom and Juan M. Vaquerizas   

*CHESS enables quantitative comparison of chromatin contact data and automatic feature extraction*    

Nature Genetics **52**, 1247-1255 (2020).


Documentation
* [CHESS Github page](https://github.com/vaquerizaslab/chess)
* [CHESS Documentation](https://chess-hic.readthedocs.io/en/latest/)
* [CHESS Issues](https://github.com/vaquerizaslab/chess/issues)


Important Notes
* Module Name: CHESS (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **CHESS\_HOME**  installation directory
	+ **CHESS\_BIN**       executable directory
	+ **CHESS\_DATA**    sample data folder



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. 
  
Sample session on a GPU node:



```

[user@biowulf ~]$ **sinteractive --mem=16g --cpus-per-task=4 --gres=lscratch:10**
[user@cn2379 ~]$ **module load chess**
[+] Loading singularity  3.8.2  on cn0890
[+] Loading chess  0.3.7

```

Generate a BED file to be used as input by the "chess sim" command:

```

[user@cn2379 ~]$ **chess pairs hg38 3000000 100000 ./hg38\_chr2\_3mb\_win\_100kb\_step\_v2.bed --chromosome chr2**
2021-10-04 15:24:32,777 INFO Running '/chess-0.3.7/envs/chess/bin/chess pairs hg38 3000000 100000 ./hg38_chr2_3mb_win_100kb_step_v2.bed --chromosome chr2'
2021-10-04 15:24:34,364 INFO CHESS version: 0.3.7
2021-10-04 15:24:34,364 INFO FAN-C version: 0.9.21
2021-10-04 15:24:34,374 INFO Finished '/chess-0.3.7/envs/chess/bin/chess pairs hg38 3000000 100000 ./hg38_chr2_3mb_win_100kb_step_v2.bed --chromosome chr2'
[user@cn2379 ~]$ **head hg38\_chr2\_3mb\_win\_100kb\_step\_v2.bed** 
chr2    0       3000000 chr2    0       3000000 0       .       +       +
chr2    100000  3100000 chr2    100000  3100000 1       .       +       +
chr2    200000  3200000 chr2    200000  3200000 2       .       +       +
chr2    300000  3300000 chr2    300000  3300000 3       .       +       +
chr2    400000  3400000 chr2    400000  3400000 4       .       +       +
chr2    500000  3500000 chr2    500000  3500000 5       .       +       +
chr2    600000  3600000 chr2    600000  3600000 6       .       +       +
chr2    700000  3700000 chr2    700000  3700000 7       .       +       +
chr2    800000  3800000 chr2    800000  3800000 8       .       +       +
chr2    900000  3900000 chr2    900000  3900000 9       .       +       +

```

Clone the chess repository, which contains sample data for running the "chess sim" command:

```
 
[user@cn2379 ~]$ **git clone https://github.com/vaquerizaslab/chess** 
[user@cn2379 ~]$ **ls chess/examples/dlbcl/\*hic**
chess/examples/dlbcl/ukm_control_fixed_le_25kb_chr2.hic
chess/examples/dlbcl/ukm_patient_fixed_le_25kb_chr2.hic

```

Run the comparison between two sets of chromatin contact data:

```

[user@cn2379 ~]$ **chess sim chess/examples/dlbcl/ukm\_control\_fixed\_le\_25kb\_chr2.hic chess/examples/dlbcl/ukm\_patient\_fixed\_le\_25kb\_chr2.hic hg38\_chr2\_3mb\_win\_100kb\_step\_v2.bed chr2\_test\_results.tsv**
2021-10-04 15:37:35,790 INFO CHESS version: 0.3.7
2021-10-04 15:37:35,791 INFO FAN-C version: 0.9.21
2021-10-04 15:37:35,792 INFO Loading reference contact data
2021-10-04 15:37:52,690 INFO Loading query contact data
2021-10-04 15:38:15,418 INFO Loading region pairs
2021-10-04 15:38:15,432 INFO Launching workers
2021-10-04 15:38:15,497 INFO Submitting pairs for comparison
2021-10-04 15:39:31,546 INFO Could not compute similarity for 105 region pairs.This can be due to faulty coordinates, too smallregion sizes or too many unmappable bins
2021-10-04 15:39:32,514 INFO Finished '/chess-0.3.7/envs/chess/bin/chess sim chess/examples/dlbcl/ukm_control_fixed_le_25kb_chr2.hic chess/examples/dlbcl/ukm_patient_fixed_le_25kb_chr2.hic hg38_chr2_3mb_win_100kb_step_v2.bed chr2_test_results.tsv'
Closing remaining open files:chess/examples/dlbcl/ukm_control_fixed_le_25kb_chr2.hic...donechess/examples/dlbcl/ukm_patient_fixed_le_25kb_chr2.hic...done

```

This command outputs the file chr2\_test\_results.tsv:

```

[user@cn2379 ~]$ **head chr2\_test\_results.tsv**
ID      SN      ssim    z_ssim
0       0.2960129270866764      0.3175164525478461      -0.5045618358016883
1       0.29579964403072184     0.32049616624170313     -0.48069735793965745
2       0.30973377588805523     0.3081302861251631      -0.5797354879768108
3       0.3216506307433863      0.27912370967940187     -0.8120486806403446
4       0.34692871691024163     0.27788389871239105     -0.8219783062403915
5       0.3524778796482955      0.26193620236586895     -0.9497031434274604
6       0.3735910245839484      0.25348665707859275     -1.017375412366634
7       0.3747548713752975      0.2497188855329951      -1.0475514325743729
8       0.3784210151741191      0.2668178225320114      -0.9106063279763916

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. chess.sh). For example:



```

#!/bin/bash
set -e
module load chess
chess sim <control_file>.hic <patient_file>.hic <chess_pairs_output>.bed results.tsv

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
**sbatch chess.sh**
```







