

document.querySelector('title').textContent = 'OrthoFinder: phylogenetic orthology inference for comparative genomics ';
**OrthoFinder: phylogenetic orthology inference for comparative genomics** 


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



OrthoFinder is an accurate and comprehensive platform for comparative genomics.
It finds orthogroups and orthologs, infers rooted gene trees for all orthogroups
and identifies all of the gene duplication events in those gene trees.



Documentation
* [FLAIR GitHub page](https://github.com/BrooksLabUCSC/flair)
* [FLAIR Manual](https://flair.readthedocs.io/en/latest/)


### References:


* David M. Emms & Steven Kelly   

 *OrthoFinder: phylogenetic orthology inference for comparative genomics*   

[Genome Biology volume 20, Article number: 238 (2019)](https://link.springer.com/article/10.1186/s13059-019-1832-y)


Important Notes
* Module Name: OrthoFinder (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **ORTHOFINDER\_HOME**  installation directory
	+ **ORTHOFINDER\_BIN**       executable directory
	+ **ORTHOFINDER\_SRC**       source code directory
	+ **ORTHOFINDER\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=4 --mem=4g --gres=lscratch:10**
[user@cn3200 ~]$ **module load orthofinder**
[+] Loading singularity  3.8.5-1  on cn0883
[+] Loading orthofinder  2.5.4
[user@cn3200 ~]$ **orthofinder -h**

OrthoFinder version 2.5.4 Copyright (C) 2014 David Emms

SIMPLE USAGE:
Run full OrthoFinder analysis on FASTA format proteomes in <dir>
  orthofinder [options] -f <dir>

Add new species in <dir1> to previous run in <dir2> and run new analysis
  orthofinder [options] -f <dir1> -b <dir2>

OPTIONS:
 -t <int>        Number of parallel sequence search threads [Default = 72]
 -a <int>        Number of parallel analysis threads
 -d              Input is DNA sequences
 -M <txt>        Method for gene tree inference. Options 'dendroblast' & 'msa'
                 [Default = dendroblast]
 -S <txt>        Sequence search program [Default = diamond]
                 Options: blast, diamond, diamond_ultra_sens, blast_gz, mmseqs, blast_nucl
 -A <txt>        MSA program, requires '-M msa' [Default = mafft]
                 Options: mafft, muscle
 -T <txt>        Tree inference method, requires '-M msa' [Default = fasttree]
                 Options: fasttree, raxml, raxml-ng, iqtree
 -s <file>       User-specified rooted species tree
 -I <int>        MCL inflation parameter [Default = 1.5]
 -x <file>       Info for outputting results in OrthoXML format
 -p <dir>        Write the temporary pickle files to <dir>
 -1              Only perform one-way sequence search
 -X              Don't add species names to sequence IDs
 -y              Split paralogous clades below root of a HOG into separate HOGs
 -z              Don't trim MSAs (columns>=90% gap, min. alignment length 500)
 -n <txt>        Name to append to the results directory
 -o <txt>        Non-default results directory
 -h              Print this help text

WORKFLOW STOPPING OPTIONS:
 -op             Stop after preparing input files for BLAST
 -og             Stop after inferring orthogroups
 -os             Stop after writing sequence files for orthogroups
                 (requires '-M msa')
 -oa             Stop after inferring alignments for orthogroups
                 (requires '-M msa')
 -ot             Stop after inferring gene trees for orthogroups

WORKFLOW RESTART COMMANDS:
 -b  <dir>         Start OrthoFinder from pre-computed BLAST results in <dir>
 -fg <dir>         Start OrthoFinder from pre-computed orthogroups in <dir>
 -ft <dir>         Start OrthoFinder from pre-computed gene trees in <dir>

LICENSE:
 Distributed under the GNU General Public License (GPLv3). See License.md

CITATION:
 When publishing work that uses OrthoFinder please cite:
 Emms D.M. & Kelly S. (2019), Genome Biology 20:238

 If you use the species tree in your work then please also cite:
 Emms D.M. & Kelly S. (2017), MBE 34(12): 3267-3278
 Emms D.M. & Kelly S. (2018), bioRxiv https://doi.org/10.1101/267914
[user@cn3200 ~]$ **git clone https://github.com/davidemms/OrthoFinder**
[user@cn3200 ~]$ **orthofinder -f ./OrthoFinder/ExampleData/** 

OrthoFinder version 2.5.4 Copyright (C) 2014 David Emms

2022-08-08 16:27:10 : Starting OrthoFinder 2.5.4
56 thread(s) for highly parallel tasks (BLAST searches etc.)
7 thread(s) for OrthoFinder algorithm

Checking required programs are installed
----------------------------------------
Test can run "mcl -h" - ok
Test can run "fastme -i /gpfs/gsfs7/users/user/OrthoFinder-2.5.4/ExampleData/OrthoFinder/Results_Aug08_8/WorkingDirectory/SimpleTest.phy -o /gpfs/gsfs7/users/user/OrthoFinder-2.5.4/ExampleData/OrthoFinder/Results_Aug08_8/WorkingDirectory/SimpleTest.tre" - ok

Dividing up work for BLAST for parallel processing
--------------------------------------------------
2022-08-08 16:27:11 : Creating diamond database 1 of 4
2022-08-08 16:27:11 : Creating diamond database 2 of 4
2022-08-08 16:27:11 : Creating diamond database 3 of 4
2022-08-08 16:27:11 : Creating diamond database 4 of 4

Running diamond all-versus-all
------------------------------
Using 56 thread(s)
2022-08-08 16:27:11 : This may take some time....
2022-08-08 16:27:18 : Done all-versus-all sequence search

Running OrthoFinder algorithm
-----------------------------
2022-08-08 16:27:18 : Initial processing of each species
2022-08-08 16:27:18 : Initial processing of species 2 complete
2022-08-08 16:27:18 : Initial processing of species 3 complete
2022-08-08 16:27:18 : Initial processing of species 0 complete
2022-08-08 16:27:18 : Initial processing of species 1 complete
2022-08-08 16:27:20 : Connected putative homologues
2022-08-08 16:27:20 : Written final scores for species 2 to graph file
2022-08-08 16:27:20 : Written final scores for species 1 to graph file
2022-08-08 16:27:20 : Written final scores for species 0 to graph file
2022-08-08 16:27:20 : Written final scores for species 3 to graph file
2022-08-08 16:27:21 : Ran MCL

Writing orthogroups to file
---------------------------
OrthoFinder assigned 2218 genes (81.2% of total) to 606 orthogroups. Fifty percent of all genes were in orthogroups with 4 or more genes (G50 was 4) and were contained in the largest 279 orthogroups (O50 was 279). There were 268 orthogroups with all species present and 245 of these consisted entirely of single-copy genes.

2022-08-08 16:27:22 : Done orthogroups

Analysing Orthogroups
=====================

Calculating gene distances
--------------------------
2022-08-08 16:27:25 : Done

Inferring gene and species trees
--------------------------------
2022-08-08 16:27:26 : Done 0 of 325
2022-08-08 16:27:26 : Done 100 of 325
2022-08-08 16:27:27 : Done 200 of 325

268 trees had all species present and will be used by STAG to infer the species tree

Best outgroup(s) for species tree
---------------------------------
2022-08-08 16:27:32 : Starting STRIDE
2022-08-08 16:27:32 : Done STRIDE
Observed 2 well-supported, non-terminal duplications. 2 support the best roots and 0 contradict them.
Best outgroups for species tree:
  Mycoplasma_hyopneumoniae
  Mycoplasma_genitalium, Mycoplasma_gallisepticum
  Mycoplasma_agalactiae

WARNING: Multiple potential species tree roots were identified, only one will be analyed.

Reconciling gene trees and species tree
---------------------------------------
Outgroup: Mycoplasma_hyopneumoniae
2022-08-08 16:27:32 : Starting Recon and orthologues
2022-08-08 16:27:32 : Starting OF Orthologues
2022-08-08 16:27:33 : Done 0 of 325
2022-08-08 16:27:33 : Done 100 of 325
2022-08-08 16:27:33 : Done 200 of 325
2022-08-08 16:27:33 : Done 300 of 325
2022-08-08 16:27:34 : Done OF Orthologues

Writing results files
=====================
2022-08-08 16:27:35 : Done orthologues

Results:
    /data/user/OrthoFinder-2.5.4/ExampleData/OrthoFinder/Results_Aug08_8/

CITATION:
 When publishing work that uses OrthoFinder please cite:
 Emms D.M. &iamp Kelly S. (2019), Genome Biology 20:238

 If you use the species tree in your work then please also cite:
 Emms D.M. &iamp Kelly S. (2017), MBE 34(12): 3267-3278
 Emms D.M. &iamp Kelly S. (2018), bioRxiv https://doi.org/10.1101/267914

```

End the interactive session:

```

[user@cn3200 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





