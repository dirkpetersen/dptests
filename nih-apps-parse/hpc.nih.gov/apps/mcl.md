

document.querySelector('title').textContent = 'Markov Clustering (MCL): a cluster algorithm for graphs ';
**Markov Clustering (MCL): a cluster algorithm for graphs** 


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



MCL implements Markov cluster algorithm. Among its applications 
is the assignment of proteins into families based on precomputed sequence 
similarity information. This approach 
does not suffer from the problems that normally hinder other protein sequence 
clustering algorithms, such as the presence of multi-domain proteins, 
promiscuous domains and fragmented proteins.



### References:


* Stijn van Dongen  

 *Graph Clustering by Flow Simulation*  

[PhD thesis,](http://www.library.uu.nl/digiarchief/dip/diss/1895620/inhoud.htm)  University of Utrecht, May 2000. 
* Enright A.J., Van Dongen S., Ouzounis C.A.  

*An efficient algorithm for large-scale detection of protein families,*   

[Nucleic Acids Research](https://academic.oup.com/nar/article/30/7/1575/2376029)**30**(7):1575-1584 (2002).


Documentation
* [MCL Home page](https://www.micans.org/mcl)
* [MCL Github page](https://github.com/micans/mcl)


Important Notes
* Module Name: MCL (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **MCL\_HOME**  installation directory
	+ **MCL\_BIN**       executable directory
	+ **MCL\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive** 
[user@cn3200 ~]$ **module load mcl** 
[+] Loading MCL 14-137  ...
[user@cn3200 ~]$ **cp $MCL\_DATA/\* .**
[user@cn3200 ~]$ **mcl cathat --abc -o out2.cathat** 
[mcl] new tab created
[mcl] pid 10660
 ite ------  chaos  time hom(avg,lo,hi) m-ie m-ex i-ex fmv
  1  ......   0.47  0.00 0.87/0.80/0.95 1.33 1.33 1.33 100
  2  ......   0.53  0.00 0.86/0.69/0.95 1.24 1.21 1.67 100
  3  ......   0.35  0.00 0.95/0.88/1.00 1.00 0.67 1.14 100
  4  ......   0.44  0.00 0.94/0.88/1.00 1.08 0.68 0.81 100
  5  ......   0.24  0.00 0.89/0.78/1.00 1.17 0.67 0.57 100
  6  ......   0.20  0.00 0.89/0.80/0.99 0.92 0.92 0.57 100
  7  ......   0.12  0.00 0.95/0.95/0.95 0.92 0.92 0.57 100
  8  ......   0.20  0.00 0.92/0.85/1.00 0.92 0.92 0.57 100
  9  ......   0.25  0.00 0.88/0.76/1.00 0.92 0.69 0.43 100
 10  ......   0.15  0.00 0.93/0.85/1.00 0.90 0.90 0.43 100
 11  ......   0.02  0.00 0.99/0.98/1.00 0.90 0.90 0.43 100
 12  ......   0.00  0.00 1.00/1.00/1.00 0.90 0.90 0.43 100
 13  ......   0.00  0.00 1.00/1.00/1.00 0.90 0.60 0.29 100
[mcl] jury pruning marks: <100,99,99>, out of 100
[mcl] jury pruning synopsis: <99.6 or perfect> (cf -scheme, -do log)
[mcl] output is in out.cathat
[mcl] 2 clusters found
[mcl] output is in out.cathat
...
...


```

End the interactive session:

```

[user@cn3200 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





