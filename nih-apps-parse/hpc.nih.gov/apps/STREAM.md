

document.querySelector('title').textContent = 'STREAM: Single-cell Trajectories Reconstruction, Exploration And Mapping of omics data ';
**STREAM: Single-cell Trajectories Reconstruction, Exploration And Mapping of omics data** 


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



STREAM is an interactive pipeline capable of disentangling and visualizing 
complex branching trajectories from both single-cell transcriptomic and
epigenomic data.



### References:


* Huidong Chen, Luca Albergante, Jonathan Y Hsu, Caleb A Lareau, Giosue Lo Bosco, Jihong Guan,
Shuigeng Zhou, Alexander N Gorban, Daniel E Bauer, Martin J Aryee, David M Langenau, Andrei
Zinovyev, Jason D Buenrostro, Guo-Cheng Yuan, Luca Pinello   

*STREAM: Single-cell Trajectories Reconstruction, Exploration And Mapping of omics data*   

[bioRxiv](https://www.biorxiv.org/content/early/2018/04/18/302554.1)  April 18, 2018, doi:10.1101/302554.


Documentation
* [STREAM on GitHub](https://github.com/pinellolab/STREAM )
* [STREAM Home Page](http://stream.pinellolab.org/)


Important Notes
* Module Name: STREAM (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **STREAM\_HOME**  STREAM installation directory
	+ **STREAM\_BIN**       STREAM executable directory
	+ **STREAM\_DATA**    STREAM sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$  **sinteractive --mem=16g -c14**
[user@@cn3107 ~]$ **module load STREAM**
[user@@cn3107 ~]$ **cp $STREAM\_DATA/\* .** 
[user@@cn3107 ~]$ **stream -m data\_Guo.tsv.gz -l cell\_label.tsv.gz -c cell\_label\_color.tsv.gz -g gene\_list.tsv.gz -s all**
   _____ _______ _____  ______          __  __
  / ____|__   __|  __ \|  ____|   /\   |  \/  |
 | (___    | |  | |__) | |__     /  \  | \  / |
  \___ \   | |  |  _  /|  __|   / /\ \ | |\/| |
  ____) |  | |  | | \ \| |____ / ____ \| |  | |
 |_____/   |_|  |_|  \_\______/_/    \_\_|  |_|


- Single-cell Trajectory Reconstruction and Mapping -

[Luca Pinello & Huidong Chen 2018, send bugs, suggestions or comments to lucapinello AT gmail DOT com]

Version 0.1.0

Loading input data...
Input: 271 cells, 175 genes
Filtering genes...
After filtering out low-expressed genes:
271 cells, 172 genes
Number of CPUs being used: 56
Reducing dimension...
Structure Learning...
Clustering...
Minimum Spanning Tree...
Number of initial branches: 5
Elastic Principal Graph...
[1] "Constructing tree 1 of 1 / Subset 1 of 1"
[1] "Computing EPG with 50 nodes on 271 points and 3 dimensions"
[1] "Using a single core"
Nodes = 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49
BARCODE ENERGY  NNODES  NEDGES  NRIBS   NSTARS  NRAYS   NRAYS2  MSE     MSEP    FVE     FVEP    UE      UR      URN     URN2    URSD
2||50   0.000668        50      49      44      2       0       0       0.0002886       0.0002736       0.974   0.9754  0.000344        3.539e-05       0.001769        0.08846 0
13.16 sec elapsed
[[1]]

Number of branches after initial ElPiGraph: 5
Collasping small branches ...
[1] "Constructing tree 1 of 1 / Subset 1 of 1"
[1] "Computing EPG with 80 nodes on 271 points and 3 dimensions"
[1] "Using a single core"
Nodes = 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79
BARCODE ENERGY  NNODES  NEDGES  NRIBS   NSTARS  NRAYS   NRAYS2  MSE     MSEP    FVE     FVEP    UE      UR      URN     URN2    URSD
2||80   0.0004009       80      79      74      2       0       0       0.0001989       0.0001921       0.9821  0.9827  0.0001738       2.82e-05        0.002256        0.1805  0
4.812 sec elapsed
Number of branches after optimization: 5
Extending leaves with additional nodes ...
Number of branches after extension: 5
Visualizing cells...
Flat tree plots...
/opt/conda/lib/python2.7/site-packages/networkx/drawing/nx_pylab.py:520: MatplotlibDeprecationWarning: The is_string_like function was deprecated in version 2.1.
  if not cb.is_string_like(edge_color) \
/opt/conda/lib/python2.7/site-packages/networkx/drawing/nx_pylab.py:524: MatplotlibDeprecationWarning: The is_string_like function was deprecated in version 2.1.
  for c in edge_color]):
/opt/conda/lib/python2.7/site-packages/networkx/drawing/nx_pylab.py:722: MatplotlibDeprecationWarning: The is_string_like function was deprecated in version 2.1.
  if not cb.is_string_like(label):
Subway map and Stream plots...
Visulizing genes...
Finished computation...


```

**NOTE:** in order to download more sample data, use the
 [Dropbox link](https://www.dropbox.com/sh/xnw9ro22bgrz2pa/AADQWmyCjUekg3hudvhsrAWka?dl=0).   

             The sample data available at the [STREAM GitHub page](https://github.com/pinellolab/STREAM ) are corrupt!!!   


