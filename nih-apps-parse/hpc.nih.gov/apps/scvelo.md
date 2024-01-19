

document.querySelector('title').textContent = 'scVelo: RNA velocity of single cell generalized through dynamical modeling.';
**scVelo: RNA velocity of single cell generalized through dynamical modeling.**


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



scVelo is a method to describe the rate of gene expression change for an individual gene 
at a given time point based on the ratio of its spliced and unspliced messenger RNA (mRNA). 
It avoids errors in the velocity estimates by solving the full transcriptional dynamics 
of splicing kinetics using a likelihood-based dynamical model. 
This generalizes RNA velocity to systems with transient cell states, 
which are common in development and in response to perturbations.



### References:


* Volker Bergen, Marius Lange, Stefan Peidli, F. Alexander Wolf and Fabian J. Theis   

*Generalizing RNA velocity to transient cell states through dynamical modeling*    

[Nature Biotechnology](https://www.nature.com/articles/s41587-020-0591-3) **38**,
1408-1414 (2020).


Documentation
* [scVelo GitHub page](https://github.com/theislab/scvelo)
* [scVelo Documentation](https://scvelo.readthedocs.io)


Important Notes
* Module Name: scVelo (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **SCVELO\_HOME**  installation directory
	+ **SCVELO\_BIN**       executable directory
	+ **SCVELO\_SRC**       source code directory
	+ **SCVELO\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=4g**
[user@cn0911 ~]$**module load scvelo** 
[+] Loading scvelo  0.2.4  on cn0911
[+] Loading singularity  3.8.5  on cn0911
[user@cn0911 ~]$**python $SCVELO\_SRC/tests/test\_basic.py**
[user@cn0911 ~]$**python** 
Python 3.8.12 (default, Oct 12 2021, 13:49:34)
[GCC 7.5.0] :: Anaconda, Inc. on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> **import numpy as np** 
>>> **import matplotlib.pyplot as pl** 
>>> **import scvelo as scv** 

>>> **scv.settings.set\_figure\_params('scvelo', dpi\_save=200, dpi=80, transparent=True)** 
>>> **scv.settings.plot\_prefix = 'scvelo\_fig2\_'** 
>>> **scv.settings.verbosity = 2** 

>>> **adata = scv.datasets.dentategyrus()** 
100%|██████████████████████████████████████████████████████████████████████████████████████| 23.7M/23.7M [00:00<00:00, 71.1MB/s]

>>> **scv.pp.filter\_and\_normalize(adata, min\_shared\_cells=20, n\_top\_genes=2000)** 
Filtered out 11835 genes that are detected in less than 20 cells (shared).
Normalized count data: X, spliced, unspliced.
Extracted 2000 highly variable genes.
Logarithmized X
>>> **scv.pp.moments(adata, n\_neighbors=30, n\_pcs=30)** 
computing neighbors
    finished (0:00:13)
computing moments based on connectivities
    finished (0:00:00)
>>> **scv.tl.velocity(adata, vkey='steady\_state\_velocity', mode='steady\_state')** 
computing velocities
    finished (0:00:00)
>>> **scv.tl.velocity\_graph(adata, vkey='steady\_state\_velocity')** 
computing velocity graph (using 1/56 cores)
    finished (0:00:04)
>>> **scv.tl.recover\_dynamics(adata)** 
recovering dynamics (using 1/56 cores)
    finished (0:06:25)
>>> **scv.tl.velocity(adata, mode='dynamical', vkey='dynamical\_velocity')** 
computing velocities
    finished (0:00:02)
>>> **scv.tl.velocity\_graph(adata, vkey='dynamical\_velocity', variance\_stabilization=True)** 
computing velocity graph (using 1/56 cores)
    finished (0:00:08)
>>> **scv.pl.velocity\_embedding\_stream(adata, vkey='dynamical\_velocity')** 
computing velocity embedding
    finished (0:00:08)

```


![](scvelo/Figure_1.png)



End the interactive session:

```

[user@cn0911 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





