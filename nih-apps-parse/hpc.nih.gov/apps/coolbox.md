

document.querySelector('title').textContent = 'CoolBox: a flexible toolkit for visual analysis of genomics data';
**CoolBox: a flexible toolkit for visual analysis of genomics data**


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



CoolBox is an open-source, user-friendly toolkit for visual analysis of genomics data. 
It is highly compatible with the Python ecosystem and customizable with a well-designed user interface. 
It can bed used, for example, to produce high-quality genome track
plots or fetch commonly used genomic data files with a Python script or command
line, to explore genomic data interactively within Jupyter environment or web browser.



### References:


* Weize Xu, Quan Zhong, Da Lin, Ya Zuo, Jinxia Dai, Guoliang Li and Gang Cao,   

*CoolBox: a flexible toolkit for visual analysis of genomics data*   

[BMC Bioinformatics (2021) 22:489](https://link.springer.com/article/10.1186/s12859-021-04408-w)


Documentation
* [CoolBox Github page](https://gangcaolab.github.io/CoolBox)
* [CoolBox Documentation](https://gangcaolab.github.io/CoolBox/index.html)


Important Notes
* Module Name: coolbox (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **CB\_HOME**  installation directory
	+ **CB\_BIN**       executable directory
	+ **CB\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive** 
[user@cig 3335 ~]$ **module load coolbox**
[+] Loading singularity  3.10.5  on cn3335
[+] Loading coolbox/0.3.8    
[user@cn3335 ~]$ **mkdir /data/$USER/coolbox &&iamp cd /data/$USER/coolbox**
[user@cn3335 ~]$ **cp $CB\_DATA/\* .**
[user@cn3335 ~]$ **wget https://github.com/GangCaoLab/CoolBox/archive/refs/tags/0.3.8.tar.gz**
[user@cn3335 ~]$ **tar -zxf 0.3.8.tar.gz &&iamp rm -f 0.3.8.tar.gz &&iamp cd CoolBox-0.3.8**
[user@cn3335 ~]$ **scripts/coolbox add XAxis - add Cool ../cool\_chr9\_4000000\_6000000.mcool - add Title "cool" - add BAMCov ../bam\_chr9\_4000000\_6000000.bam - add Title "bam" - add BigWig ../bigwig\_chr9\_4000000\_6000000.bw - goto "chr9:4000000-6000000" - plot test\_coolbox.png**

```

![](coolbox/test_coolbox.png)




