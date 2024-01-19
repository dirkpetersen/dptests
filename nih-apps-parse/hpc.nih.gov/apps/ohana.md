

document.querySelector('title').textContent = 'Ohana on Biowulf';
Ohana on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
 |



Ohana is a suite of software for analyzing population structure and admixture history using unsupervised learning methods. We construct statistical models to infer individual clustering from which we identify outliers for selection analyses.



Documentation
* [Ohana Main Site](https://github.com/jade-cheng/ohana)


Important Notes
* Module Name: ohana (see [the modules page](/apps/modules.html) for more information)
 * singlethreaded



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load ohana**
[user@cn3144 ~]$ **convert ped2dgm sample.ped g.dgm**
[user@cn3144 ~]$ **qpas g.dgm -k 4 -qo q.matrix -fo f.matrix -mi 5**
[user@cn3144 ~]$ **python $OHANA\_HOME/tools/plot\_q.py q.matrix q-bar-chart.pdf**
[user@cn3144 ~]$ **convert cov2nwk f.matrix tree.nwk**
[user@cn3144 ~]$ **convert nwk2svg tree.nwk tree.svg**
[user@cn3144 ~]$ **selscan g.dgm f.matrix c.matrix > lle-ratios.txt**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





