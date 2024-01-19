

document.querySelector('title').textContent = 'LTSOFT: analysis of case–control association studies with known risk variants';
**LTSOFT: analysis of case–control association studies with known risk variants**


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



The LTSOFT application implements a new approach 
to using information from known associated variants 
when conducting disease association studies.
The approach is based in part on the classical technique 
of liability threshold modeling and performs estimation
of model parameters for each known variant 
while accounting for the published disease prevalence from
the epidemiological literature.



### References:


* Noah Zaitlen, Bogdan Pasaniuc, Nick Patterson, et al.,   

*Analysis of case–control association studies with known risk variants*  

[Bioinformatics](https://academic.oup.com/bioinformatics/article/28/13/1729/235461?login=true)Volume 28, Issue 13, July 2012, Pages 1729–1737, https://doi.org/10.1093/bioinformatics/bts259.
* Noah Zaitlen, Sara Lindstrom, Bogdan Pasaniuc, et al.,  

*Informed Conditioning on Clinical Covariates Increases Power in Case-Control Association Studies*  

[PLOS Genetics](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003032)November 2012 | Volume 8 | Issue 11 | e1003032.


Documentation
* [LTSOFT Home page](https://www.hsph.harvard.edu/alkes-price/software/)


Important Notes
* Module Name: ltsoft (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **LTSOFT\_HOME**  installation directory
	+ **LTSOFT\_BIN**       executable directory
	+ **LTSOFT\_SRC**       source code directory
	+ **LTSOFT\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive** 
[user@cig 3335 ~]$ **module load ltsoft**
[+] Loading singularity  3.10.5  on cn4326
[+] Loading ltsoft  4.0
[user@cn4326 ~]$ **ls $LTSOFT\_BIN**
convertf  gcta64  ltfit  LTFIT  LTMLM  ltprev  ltpub  ltscore  LTSCORE  ltscore.perl 
[user@cn4326 ~]$ **ltfit**
usage:  ltfit  < parameter-file >
[user@cn4326 ~]$ **ltprev**
usage:  ltprev < parameter-file >
[user@cn4326 ~]$ **ltpub**
usage: ltpub < in-file > < outfile >
[user@cn4326 ~]$ **ltscore**
usage:  ltscore 
[user@cn4326 ~]$ **ltscore.perl --help**
...

Usage: ltscore.perl [-OPTIONS [-MORE\_OPTIONS]] [--] [PROGRAM\_ARG1 ...]

The following single-character options are accepted:
 With arguments: -g -i -n -c -o -s -p

Options may be merged together. -- stops processing of options.
Space is not required between options and their arguments.
...

```





