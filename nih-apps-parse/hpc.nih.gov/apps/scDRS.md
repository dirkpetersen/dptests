

document.querySelector('title').textContent = 'scDRS: disease relevance scores for individual cells in single-cell RNA-seq data';
**scDRS: disease relevance scores for individual cells in single-cell RNA-seq data**


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



The scDRS application implements an approach that links scRNA-seq 
with polygenic disease risk at single-cell resolution, independent of annotated cell types. 
scDRS identifies cells exhibiting excess expression across disease-associated genes 
implicated by genome-wide association studies (GWASs). 
Genes whose expression was correlated with the scDRS score across cells 
(reflecting coexpression with GWAS disease-associated genes) 
were strongly enriched for gold-standard drug target and Mendelian disease genes.



### References:


* Martin Jinye Zhang, Kangcheng Hou, Kushal K. Dey, Saori Sakaue, Karthik A. Jagadeesh, Kathryn Weinand, 
Aris Taychameekiatchai, Poorvi Rao, Angela Oliveira Pisco, James Zou, Bruce Wang, Michael Gandal, 
Soumya Raychaudhuri, Bogdan Pasaniuc &iamp Alkes L. Price   

*Polygenic enrichment distinguishes disease associations of individual cells in single-cell RNA-seq data*    

[Nature Genetics](https://www.nature.com/articles/s41588-022-01167-z)  **54**, pages 1572–1580 (2022).


Documentation
* [scDRS GitHub page](https://github.com/martinjzhang/scDRS)
* [scDRS Documentation](https://martinjzhang.github.io/scDRS/notebooks/quickstart.html)


Important Notes
* Module Name: scdrs (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **SCDRS\_HOME**  installation directory
	+ **SCDRS\_BIN**       executable directory
	+ **SCDRS\_SRC**       source code directory
	+ **SCDRS\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=4g**
[user@cn0911 ~]$**module load scdrs** 
[+] Loading singularity  3.10.5  on cn4185
[+] Loading scdrs  1.02 
[user@cn0911 ~]$**mkdir /data/$USER/scDRS && cd /data/$USER/scDRS**
[user@cn0911 ~]$**wget https://figshare.com/ndownloader/files/34300925 -O data.zip**
[user@cn0911 ~]$**unzip data.zip && \
 mkdir -p data/ && \
 mv single\_cell\_data/zeisel\_2015/\* data/ && \
 rm data.zip && rm -r single\_cell\_data**
[user@cn0911 ~]$**scdrs compute-score \
 --h5ad-file data/expr.h5ad \
 --h5ad-species mouse \
 --gs-file data/geneset.gs \
 --gs-species mouse \
 --cov-file data/cov.tsv \
 --flag-filter-data True \
 --flag-raw-count True \
 --flag-return-ctrl-raw-score False \
 --flag-return-ctrl-norm-score True \
 --out-folder data/**

Loading data:
--h5ad-file loaded: n_cell=3005, n_gene=13572 (sys_time=1.7s)
First 3 cells: ['1772071015_C02', '1772071017_G12', '1772071017_A05']
First 5 genes: ['Tspan12', 'Tshz1', 'Fnbp1l', 'Adamts15', 'Cldn12']
--cov-file loaded: covariates=['const', 'n_genes'] (sys_time=1.7s)
First 5 values for 'const': [1, 1, 1, 1, 1]
First 5 values for 'n_genes': [4848, 4685, 6028, 5824, 4701]
--gs-file loaded: n_trait=29 (sys_time=1.8s)
Print info for first 3 traits:
First 3 elements for 'PASS_ADHD_Demontis2018': ['St3gal3', 'Kdm4a', 'Ptprf'], [6.4588, 6.2164, 5.8681]
First 3 elements for 'PASS_Alzheimers_Jansen2019': ['Ms4a6d', 'Clu', 'Picalm'], [7.1313, 7.066, 6.6287]
First 3 elements for 'PASS_BIP_Mullins2021': ['Trank1', 'Myrf', 'Fads1'], [7.2182, 7.0688, 6.9447]

Preprocessing:

Computing scDRS score:
Computing control scores: 100%|█████████████████████████████████████████████| 1000/1000 [00:53<00:00, 18.74it/s]
Trait=PASS_ADHD_Demontis2018, n_gene=688: 2/3005 FDR<0.1 cells, 4/3005 FDR<0.2 cells (sys_time=82.3s)
Computing control scores: 100%|███████████████████████████████████████████████████| 1000/1000 [00:49<00:00, 20.27it/s]
Trait=PASS_Alzheimers_Jansen2019, n_gene=630: 0/3005 FDR<0.1 cells, 0/3005 FDR<0.2 cells (sys_time=155.8s)
Computing control scores: 100%|███████████████████████████████████████████████████| 1000/1000 [00:58<00:00, 17.24it/s]
Trait=PASS_BIP_Mullins2021, n_gene=776: 171/3005 FDR<0.1 cells, 314/3005 FDR<0.2 cells (sys_time=238.8s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:54<00:00, 18.34it/s]
Trait=PASS_BIP_Stahl2019, n_gene=716: 77/3005 FDR<0.1 cells, 186/3005 FDR<0.2 cells (sys_time=317.8s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:52<00:00, 18.89it/s]
Trait=PASS_DrinksPerWeek_Liu2019, n_gene=694: 0/3005 FDR<0.1 cells, 0/3005 FDR<0.2 cells (sys_time=394.5s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:53<00:00, 18.66it/s]
Trait=PASS_GeneralRiskTolerance_KarlssonLinner2019, n_gene=700: 173/3005 FDR<0.1 cells, 314/3005 FDR<0.2 cells (sys_time=472.5s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:54<00:00, 18.43it/s]
Trait=PASS_Insomnia_Jansen2019, n_gene=710: 0/3005 FDR<0.1 cells, 1/3005 FDR<0.2 cells (sys_time=550.9s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:55<00:00, 17.93it/s]
Trait=PASS_Intelligence_SavageJansen2018, n_gene=742: 170/3005 FDR<0.1 cells, 299/3005 FDR<0.2 cells (sys_time=631.0s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:54<00:00, 18.34it/s]
Trait=PASS_MDD_Howard2019, n_gene=716: 281/3005 FDR<0.1 cells, 434/3005 FDR<0.2 cells (sys_time=709.5s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:55<00:00, 18.07it/s]
Trait=PASS_ReactionTime_Davies2018, n_gene=733: 73/3005 FDR<0.1 cells, 191/3005 FDR<0.2 cells (sys_time=789.6s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:53<00:00, 18.81it/s]
Trait=PASS_SWB, n_gene=690: 292/3005 FDR<0.1 cells, 466/3005 FDR<0.2 cells (sys_time=867.2s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:56<00:00, 17.63it/s]
Trait=PASS_Schizophrenia_Pardinas2018, n_gene=756: 178/3005 FDR<0.1 cells, 320/3005 FDR<0.2 cells (sys_time=948.2s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:54<00:00, 18.37it/s]
Trait=PASS_SleepDuration_Dashti2019, n_gene=718: 104/3005 FDR<0.1 cells, 208/3005 FDR<0.2 cells (sys_time=1026.9s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:56<00:00, 17.66it/s]
Trait=PASS_VerbalNumericReasoning_Davies2018, n_gene=749: 157/3005 FDR<0.1 cells, 305/3005 FDR<0.2 cells (sys_time=1108.0s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:55<00:00, 18.05it/s]
Trait=PASS_Worry_Nagel2018, n_gene=733: 179/3005 FDR<0.1 cells, 285/3005 FDR<0.2 cells (sys_time=1187.8s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:56<00:00, 17.72it/s]
Trait=UKB_460K.body_BMIz, n_gene=751: 182/3005 FDR<0.1 cells, 339/3005 FDR<0.2 cells (sys_time=1269.0s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:52<00:00, 19.18it/s]
Trait=UKB_460K.body_HEIGHTz, n_gene=671: 0/3005 FDR<0.1 cells, 6/3005 FDR<0.2 cells (sys_time=1345.0s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:55<00:00, 17.86it/s]
Trait=UKB_460K.cov_EDU_COLLEGE, n_gene=743: 197/3005 FDR<0.1 cells, 338/3005 FDR<0.2 cells (sys_time=1425.2s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:54<00:00, 18.24it/s]
Trait=UKB_460K.cov_EDU_YEARS, n_gene=729: 239/3005 FDR<0.1 cells, 403/3005 FDR<0.2 cells (sys_time=1504.2s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:54<00:00, 18.27it/s]
Trait=UKB_460K.cov_SMOKING_STATUS, n_gene=725: 194/3005 FDR<0.1 cells, 360/3005 FDR<0.2 cells (sys_time=1583.4s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:55<00:00, 17.87it/s]
Trait=UKB_460K.mental_NEUROTICISM, n_gene=751: 229/3005 FDR<0.1 cells, 401/3005 FDR<0.2 cells (sys_time=1663.7s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:53<00:00, 18.65it/s]
Trait=UKB_460K.other_MORNINGPERSON, n_gene=712: 172/3005 FDR<0.1 cells, 320/3005 FDR<0.2 cells (sys_time=1741.3s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:52<00:00, 18.99it/s]
Trait=UKB_460K.repro_NumberChildrenEverBorn_Pooled, n_gene=695: 43/3005 FDR<0.1 cells, 86/3005 FDR<0.2 cells (sys_time=1817.8s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:24<00:00, 41.63it/s]
Trait=spatial_ventral, n_gene=228: 672/3005 FDR<0.1 cells, 770/3005 FDR<0.2 cells (sys_time=1859.7s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:19<00:00, 51.75it/s]
Trait=spatial_dorsal, n_gene=155: 549/3005 FDR<0.1 cells, 641/3005 FDR<0.2 cells (sys_time=1895.8s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:13<00:00, 76.27it/s]
Trait=spatial_distal, n_gene=41: 3/3005 FDR<0.1 cells, 10/3005 FDR<0.2 cells (sys_time=1923.2s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:16<00:00, 61.84it/s]
Trait=spatial_proximal, n_gene=98: 219/3005 FDR<0.1 cells, 356/3005 FDR<0.2 cells (sys_time=1954.9s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:17<00:00, 56.28it/s]
Trait=spatial_deep, n_gene=122: 484/3005 FDR<0.1 cells, 576/3005 FDR<0.2 cells (sys_time=1988.6s)
Computing control scores: 100%|███████████████████████████████████████████████████████████████| 1000/1000 [00:13<00:00, 72.87it/s]
Trait=spatial_superficial, n_gene=55: 190/3005 FDR<0.1 cells, 322/3005 FDR<0.2 cells (sys_time=2016.9s)
[user@cn0911 ~]$**exit**

```
 
End the interactive session:

```

[user@cn0911 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





