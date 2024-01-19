

document.querySelector('title').textContent = 'BioNetGen: specification and simulation of rule-based models of biochemical systems';
Lifelines: a pure Python implementation of the best parts of survival analysis.


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



Lifelines is a complete survival analysis library, written in pure Python. 
It has the benefits of easy installation, internal plotting methods, simple and intuitive API
handles right, left and interval censored data. It also contains the most popular parametric, 
semi-parametric and non-parametric models.
  

  




### References:


* Cameron Davidson-Pilon
*lifelines: survival analysis in Python*    

 [Journal of Open Source Software, 4(40), 1317. https://doi.org/10.21105/joss.
01317
1](chrome-extension://bdfcnmeidppjeaggnmidamkiddifkdib/viewer.html?file=https://joss.theoj.org/papers/10.21105/joss.01317.pdf)


Documentation
* [Lifelines Github page](https://github.com/CamDavidsonPilon/lifelines)
* [Lifelines documentation](https://lifelines.readthedocs.io/en/latest/)
* [Lifelines Video Tutorial](https://medium.com/analytics-vidhya/survival-analysis-using-lifelines-in-python-bf5eb0435dec)


Important Notes
* Module Name: lifelines (see [the modules page](/apps/modules.html) for more information)
  
* Unusual environment variables set
	+ **LIFELINES\_HOME**  installation directory
	+ **LIFELINES\_BIN**  executable directory
	+ **LIFELINES\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
[user@cn3144 ~]$ **module load lifelines**
[+] Loading singularity  3.10.5  on cn4338
[+] Loading lifelines  0.27.4

```

  

The basdic usage of the lifelines application is as follows:

```

[user@cn3144 ~]$ **sinteractive**
[user@cn3144 ~]$ **python**
>>> import lifelines
>>> from lifelines.datasets import load_waltons
>>> print(df.head())
      T  E    group
0   6.0  1  miR-137
1  13.0  1  miR-137
2  13.0  1  miR-137
3  13.0  1  miR-137
4  19.0  1  miR-137
>>> T = df['T']
>>> E = df['E']
>>> kmf = KaplanMeierFitter()
>>> kmf.fit(T, event_observed=E)

>>> kmf.survival\_function\_
 KM\_estimate
timeline
0.0 1.000000
6.0 0.993865
7.0 0.987730
9.0 0.969210
13.0 0.950690
15.0 0.938344
17.0 0.932170
19.0 0.913650
22.0 0.888957
26.0 0.858090
29.0 0.827224
32.0 0.821051
33.0 0.802531
36.0 0.790184
38.0 0.777837
41.0 0.734624
43.0 0.728451
45.0 0.672891
47.0 0.666661
48.0 0.616817
51.0 0.598125
53.0 0.554512
54.0 0.542051
56.0 0.429903
58.0 0.404981
60.0 0.311524
61.0 0.254305
62.0 0.240921
63.0 0.180690
66.0 0.160614
68.0 0.100384
69.0 0.014341
75.0 0.000000
>>> kmf.cumulative\_density\_
 KM\_estimate
timeline
0.0 0.000000
6.0 0.006135
7.0 0.012270
9.0 0.030790
13.0 0.049310
15.0 0.061656
17.0 0.067830
19.0 0.086350
22.0 0.111043
26.0 0.141910
29.0 0.172776
32.0 0.178949
33.0 0.197469
36.0 0.209816
38.0 0.222163
41.0 0.265376
43.0 0.271549
45.0 0.327109
47.0 0.333339
48.0 0.383183
51.0 0.401875
53.0 0.445488
54.0 0.457949
56.0 0.570097
58.0 0.595019
60.0 0.688476
61.0 0.745695
62.0 0.759079
63.0 0.819310
66.0 0.839386
68.0 0.899616
69.0 0.985659
75.0 1.000000

```

  

End the interactive session:

```

[user@cn3111 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





