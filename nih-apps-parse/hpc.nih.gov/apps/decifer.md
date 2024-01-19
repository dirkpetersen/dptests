

document.querySelector('title').textContent = 'DeCiFer: measuring tumor heterogenicity through clustering SNVs by their corresponding descendant cell fractions. ';
**DeCiFer: measuring tumor heterogenicity through clustering SNVs by their corresponding descendant cell fractions.** 


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



DeCiFer is an algorithm that simultaneously selects mutation multiplicities and clusters somatic single-nucleotide variants (SNVs) by their corresponding descendant cell fractions (DCF), a statistic that quantifies the proportion of cells which acquired the SNV or whose ancestors acquired the SNV. DCF is related to the commonly used cancer cell fraction (CCF) but further accounts for SNVs which are lost due to deleterious somatic copy-number aberrations (CNAs), identifying clusters of SNVs which occur in the same phylogenetic branch of tumour evolution.



### References:


* Gryte Satas, Simone Zaccaria, Mohammed El-Kebir, and Benjamin J. Raphael1   

*DeCiFering the elusive cancer cell fraction in tumor heterogeneity and evolution*   

[Cell Systems, 2021, v 12, 1004–1018](https://www.sciencedirect.com/science/article/pii/S2405471221002842).


Documentation
* [DeCiFer Github page](https://github.com/raphael-group/decifer)


Important Notes
* Module Name: decifer (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **DECIFER\_HOME**  installation directory
	+ **DECIFER\_BIN**       executable directory
	+ **DECIFER\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=24g --cpus-per-task=24 --gres=lscratch:20**
[user@cn3335 ~]$**module load decifer** 
[+] Loading singularity  3.10.3  on cn0802
[+] Loading decifer 2.1.3  ...

```

Set up references:

```

[user@biowulf]$ **decifer -h**
usage: decifer [-h] -p PURITYFILE [--betabinomial] [-i SNPFILE] [-s SEGFILE] [-v SENSITIVITY]
               [-R RESTARTS_BB] [-x SKIP] [--ccf] [-k MINK] [-K MAXK] [-r RESTARTS] [-t MAXIT]
               [-e ELBOW] [--binarysearch] [--record] [-j JOBS] [-o OUTPUT]
               [--statetrees STATETREES] [--seed SEED] [--debug] [--printallk] [--conservativeCIs]
               [--vafdevfilter VAFDEVFILTER] [--silhouette]
               INPUT

DeCiFer.

positional arguments:
  INPUT                 Input file in DeCiFer format.

optional arguments:
  -h, --help            show this help message and exit
  -p PURITYFILE, --purityfile PURITYFILE
                        File with purity of each sample (TSV file in two columns`SAMPLE PURITY`)
  --betabinomial        Use betabinomial likelihood to cluster mutations (default: binomial)
  -i SNPFILE, --snpfile SNPFILE
                        File with precisions for betabinomial fit (default: binomial likelihood)
  -s SEGFILE, --segfile SEGFILE
                        File with precisions for betabinomial fit (default: binomial likelihood)
  -v SENSITIVITY, --sensitivity SENSITIVITY
                        Sensitivity E to exclude SNPs with 0.5 - E <= BAF < 0.5, for estimating
                        betabinomial parameters (default: 0.1)
  -R RESTARTS_BB, --restarts_bb RESTARTS_BB
                        Maximum size of brute-force search, when fitting betabinomial parameters
                        (default: 1e4)
  -x SKIP, --skip SKIP  Numbers to skip in the brute-force search, when fitting betabinomial
                        parameters (default: 10)
  --ccf                 Run with CCF instead of DCF (default: False)
  -k MINK, --mink MINK  Minimum number of clusters, which must be at least 2 (default: 2)
  -K MAXK, --maxk MAXK  Maximum number of clusters (default: 12)
  -r RESTARTS, --restarts RESTARTS
                        Number of restarts (default: 20)
  -t MAXIT, --maxit MAXIT
                        Maximum number of iterations per restart (default: 200)
  -e ELBOW, --elbow ELBOW
                        Elbow sensitivity, lower values increase sensitivity (default: 0.06)
  --binarysearch        Use binary-search model selection (default: False, iterative is used; use
                        binary search when considering large numbers of clusters
  --record              Record objectives (default: False)
  -j JOBS, --jobs JOBS  Number of parallele jobs to use (default: equal to number of available
                        processors)
  -o OUTPUT, --output OUTPUT
                        Output prefix (default: ./decifer)
  --statetrees STATETREES
                        Filename of state-trees file (default: use state_trees.txt in the package)
  --seed SEED           Random-generator seed (default: None)
  --debug               single-threaded mode for development/debugging
  --printallk           Print all results for each value of K explored by DeCiFer
  --conservativeCIs     Beta: compute CIs using DCF point values assigned to cluster instead of
                        cluster likelihood function
  --vafdevfilter VAFDEVFILTER
                        Filter poorly fit SNVs with VAFs that are more than this number of
                        standard deviations away from the cluster center VAF (default 1.5)
  --silhouette          Beta: select the number of clusters using a silhouette score
...

```

Downloading mutation input file:

```

[user@cn3335 ~]$ **mkdir -p data** 
[user@cn3335 ~]$ **curl -L 'https://raw.githubusercontent.com/raphael-group/decifer-data/main/input/prostate/mutations/A12.decifer.input.tsv' > data/mutations.tsv** 

```

Downloading purity input file:

```

[user@cn3335 ~]$  **curl -L 'https://raw.githubusercontent.com/raphael-group/decifer-data/main/input/prostate/purity/A12.purity.txt' > data/purity.tsv**

```

Running DeCiFer:

```

[user@cn3335 ~]$ **decifer data/mutations.tsv -p data/purity.tsv -k 5 -K 8 -r 20 --seed 17 -j 24**
Arguments:
        input : data/mutations.tsv
        mink : 5
        maxk : 8
        maxit : 200
        purity : data/purity.tsv
        restarts : 20
        elbow : 0.06
        iterative : True
        record : False
        J : 128
        output : ./decifer
        ccf : False
        betabinomial : False
        snpfile : None
        segfile : None
        restarts_bb : 10000
        threshold : 0.1
        skip : 10
        statetrees : /opt/conda/lib/python3.9/site-packages/decifer/state_trees.txt
        debug : False
        printallk : False
        conservativeCIs : False
        vafdevfilter : 1.5
        silhouette : False
0 D
1 C
2 A
Using iterative model selection
Progress: |------------------------------| 1.2% Complete [[2022-Nov-22 16:37:38]Completed 2 for k=5 
Progress: |------------------------------| 2.5% Complete [[2022-Nov-22 16:37:38]Completed 1 for k=5
Progress: |█-----------------------------| 3.8% Complete [[2022-Nov-22 16:37:40]Completed 5 for k=5
 Progress: |█-----------------------------| 5.0% Complete [[2022-Nov-22 16:37:40]Completed 4 for k=5 
 Progress: |█-----------------------------| 6.2% Complete [[2022-Nov-22 16:37:40]Completed 3 for k=5
 Progress: |██----------------------------| 7.5% Complete [[2022-Nov-22 16:37:40]Completed 0 for k=5
 Progress: |██----------------------------| 8.8% Complete [[2022-Nov-22 16:37:48]Completed 1 for k=6
 Progress: |███---------------------------| 10.0% Complete [[2022-Nov-22 16:37:59]Completed 0 for k=7
Progress: |███---------------------------| 11.2% Complete [[2022-Nov-22 16:37:59]Completed 4 for k=7
Progress: |███---------------------------| 12.5% Complete [[2022-Nov-22 16:37:59]Completed 1 for k=7
Progress: |████--------------------------| 13.8% Complete [[2022-Nov-22 16:38:04]Completed 3 for k=6
Progress: |████--------------------------| 15.0% Complete [[2022-Nov-22 16:38:04]Completed 0 for k=6
Progress: |████--------------------------| 16.2% Complete [[2022-Nov-22 16:38:05]Completed 4 for k=6
Progress: |█████-------------------------| 17.5% Complete [[2022-Nov-22 16:38:05]Completed 2 for k=6
Progress: |█████-------------------------| 18.8% Complete [[2022-Nov-22 16:38:08]Completed 6 for k=5
...
Progress: |████████████████████████████--| 95.0% Complete [[2022-Nov-22 16:41:16]Completed 17 for k=7
Progress: |████████████████████████████--| 96.2% Complete [[2022-Nov-22 16:41:19]Completed 11 for k=6
Progress: |█████████████████████████████-| 97.5% Complete [[2022-Nov-22 16:41:33]Completed 9 for k=8
Progress: |█████████████████████████████-| 98.8% Complete [[2022-Nov-22 16:41:40]Completed 18 for k=6
Progress: |██████████████████████████████| 100.0% Complete [[2022-Nov-22 16:42:01]Completed 19 for k=8 
[Iterations: 17]]
[user@cn3335 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





