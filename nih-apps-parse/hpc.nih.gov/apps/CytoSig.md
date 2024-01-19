

document.querySelector('title').textContent = 'CytoSig: prediction of cytokine signaling activity';
**CytoSig: prediction of cytokine signaling activity**


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



CytoSig is a data-driven infrastructure hosted by the National Cancer Institute. 
CytoSig includes both a database of target genes modulated by cytokines and 
a predictive model of cytokine signaling activity and regulatory cascade from transcriptomic profiles.



### References:


* Peng Jiang, Yu Zhang, Beibei Ru, Yuan Yang, Trang Vu, Rohit Pau, Amer Mirza, Grégoire Altan-Bonnet, Lingrui Liu, Eytan Ruppin, Lalage Wakefield and Kai W. Wucherpfennig,   

*Systematic investigation of cytokine signaling activity at the tissue and single-cell levels*   

[Nature Methods](https://www.nature.com/articles/s41592-021-01274-5) 2021, VOL. **18** 1181–1191.


Documentation
* [CytoSig Home Page](https://cytosig.ccr.cancer.gov/)
* [CytoSig Github page](https://github.com/data2intelligence/CytoSig)


Important Notes
* Module Name: CytoSig (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **CYTOSIG\_HOME**  installation directory
	+ **CYTOSIG\_BIN**       executable directory
	+ **CYTOSIG\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=8g -c8 --gres=lscratch:10**
[user@cig 3335 ~]$ **module load CytoSig**
[+] Loading singularity  3.10.5  on cn3335
[+] Loading CytoSig  0.1
[user@cn3335 ~]$ **CytoSig\_run.py -h**
Usage:
CytoSig_run.py -i <input profiles> -o <output prefix> -r <randomization count, default: 1000> -a <penalty alpha, default: 10000> -e <generate excel report: 0|1, default: 0> -s <use an expanded response signature: 0|1, default: 0>


```

Download the CytoSig repository with test data to the current folder:

```

[user@cn3335 ~]$ **wget https://github.com/data2intelligence/CytoSig/archive/refs/tags/v0.1.tar.gz**
[user@cn3335 ~]$ **tar -zxf v0.1.tar.gz & rm -f v0.1.tar.gz**
[user@cn3335 ~]$ **cd CytoSig-0.1**

```

Run the unit test:

```

[user@cn3335 ~]$ **python3 -m unittest tests.prediction**
Output:
 Use permutation test with nrand = 1000
0%
10%
20%
30%
40%
50%
60%
70%
80%
90%
Use permutation test with nrand = 1000
0%
10%
20%
30%
40%
50%
60%
70%
80%
90%
.
----------------------------------------------------------------------
Ran 1 test in 3.095s

OK

```

Run SitoSig on test data:

```

[user@cn3335 ~]$ **CytoSig\_run.py -i CytoSig-0.1/tests/GSE147507.diff.gz -o output** 
Use permutation test with nrand = 1000
0%
10%
20%
30%
40%
50%
60%
70%
80%
90%

```





