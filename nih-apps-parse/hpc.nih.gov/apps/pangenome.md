

pangenome on Biowulf


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



From the pangenome user manual:




>  
>  nf-core/pangenome is a bioinformatics best-practise analysis pipeline for the rendering of a collection of sequences into a pangenome graph. Its goal is to build a graph that is locally directed and acyclic while preserving large-scale variation. Maintaining local linearity is important for interpretation, visualization, mapping, comparative genomics, and reuse of pangenome graphs.
> 
> 


Documentation
* [Manual](https://github.com/nf-core/pangenome)


Important Notes
* Module Name: pangenome (see [the modules page](/apps/modules.html) for more information)
* pangenome is a nextflow pipeline, and it can also run individually with the vg, pggb, odgi and wfmash command as shown below.

```

    [user@nc3144]$ **module load pangenome**
    [user@nc3144]$ **vg**
    
```
* Example files in `$PANGENOME_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:


Allocate an interactive session with [sinteractive](/docs/userguide.html#int)
and use as shown below. In this case we will use test data that is 
an artificial mixture of 1M human exome reads and 1M environmental metagenomic 
reads. The 50% human reads is treated as an artificial contamination and removed:



```

[user@biowulf]$ **sinteractive --mem=36g --cpus-per-task=12 --gres=lscratch:10**
salloc.exe: Pending job allocation 33247354
salloc.exe: job 33247354 queued and waiting for resources
salloc.exe: job 33247354 has been allocated resources
salloc.exe: Granted job allocation 33247354
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job
srun: error: x11: no local DISPLAY defined, skipping

[user@nc3144]$ **module load pangenome**
[user@nc3144]$ **cd /lscratch/${SLURM\_JOB\_ID}**
[user@nc3144]$ **vg**
vg: variation graph tool, version v1.40.0 "Suardi"

usage: /usr/local/bin/vg  [options]

main mapping and calling pipeline:
  -- autoindex     mapping tool-oriented index construction from interchange formats
  -- construct     graph construction
  -- rna           construct splicing graphs and pantranscriptomes
  -- index         index graphs or alignments for random access or mapping
  -- map           MEM-based read alignment
  -- giraffe       fast haplotype-aware short read alignment
  -- mpmap         splice-aware multipath alignment of short reads
  -- augment       augment a graph from an alignment
  -- pack          convert alignments to a compact coverage index
  -- call          call or genotype VCF variants
  -- help          show all subcommands

For more commands, type `vg help`.
For technical support, please visit: https://www.biostars.org/t/vg/

[user@nc3144]$ **pggb**
[user@nc3144]$ **odgi**
odgi: optimized dynamic genome/graph implementation, version v0.8.3-0-g34f006f3

usage: /usr/local/bin/odgi  [options]

Overview of available commands:
  -- bin           Binning of pangenome sequence and path information in the graph.
  -- break         Break cycles in the graph and drop its paths.
  -- build         Construct a dynamic succinct variation graph in ODGI format from a GFAv1.
  -- chop          Divide nodes into smaller pieces preserving node topology and order.
  -- cover         Cover the graph with paths.
  -- crush         Crush runs of N.
  -- degree        Describe the graph in terms of node degree.
  -- depth         Find the depth of a graph as defined by query criteria.
  -- draw          Draw previously-determined 2D layouts of the graph with diverse annotations.
  -- explode       Breaks a graph into connected components storing each component in its own file.
  -- extract       Extract subgraphs or parts of a graph defined by query criteria.
  -- flatten       Generate linearizations of a graph.
  -- flip          Flip path orientations to match the graph.
  -- groom         Harmonize node orientations.
  -- heaps         Path pangenome coverage permutations.
  -- inject        Inject BED annotations as paths.
  -- kmers         Display and characterize the kmer space of a graph.
  -- layout        Establish 2D layouts of the graph using path-guided stochastic gradient descent.
  -- matrix        Write the graph topology in sparse matrix format.
  -- normalize     Compact unitigs and simplify redundant furcations.
  -- overlap       Find the paths touched by given input paths.
  -- panpos        Get the pangenome position of a given path and nucleotide position (1-based).
  -- pathindex     Create a path index for a given graph.
  -- paths         Interrogate the embedded paths of a graph.
  -- pav           Presence/absence variants (PAVs).
  -- position      Find, translate, and liftover graph and path positions between graphs.
  -- priv          Differentially private sampling of graph subpaths.
  -- procbed       Procrustes-BED: adjust BED to match subpaths in graph.
  -- prune         Remove parts of the graph.
  -- server        Start a basic HTTP server to lift coordinates between path and pangenomic positions.
  -- similarity    Provides a sparse similarity matrix for paths or groups of paths.
  -- sort          Apply different kind of sorting algorithms to a graph.
  -- squeeze       Squeezes multiple graphs in ODGI format into the same file in ODGI format.
  -- stats         Metrics describing a variation graph and its path relationship.
  -- stepindex     Generate a step index and access the position of each step of each path once.
  -- tension       evaluate the tension of a graph helping to locate structural variants and abnormalities
  -- tips          Identifying break point positions relative to given references.
  -- unchop        Merge unitigs into a single node preserving the node order.
  -- unitig        Output unitigs of the graph.
  -- untangle      Project paths into reference-relative, to decompose paralogy relationships.
  -- validate      Validate a graph checking if the paths are consistent with the graph topology.
  -- version       Print the version of ODGI to stdout.
  -- view          Project a graph into other formats.
  -- viz           Visualize a variation graph in 1D.

[user@nc3144]$ # copy test data
[user@nc3144]$ **cp $PANGENOME\_TEST\_DATA/\* .**
[user@nc3144]$ **cp /usr/local/apps/nextflow/nextflow.config .**
[user@nc3144]$ **# run pangenome on the test data with biowulflocal profile**
[user@nc3144]$ **nextflow run nf-core/pangenome \
-r dev -profile biowulflocal \
--input DRB1-3123.fa.gz \
--n\_haplotypes 12 \
--outdir testout \
-with-singularity $PANGENOME\_PATH/pangenome.img**
[user@nc3144]$ **exit**
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. pangenome.sh) similar to the following example:



```

#! /bin/bash

module load pangenome || exit 1

if [[ ! -e DRB1-3123.fa.gz ]]; then
    cp $PANGENOME_TEST_DATA/* .
fi
cp /usr/local/apps/nextflow/nextflow.config .

nextflow run nf-core/pangenome \
-r 1.0.0 -profile biowulf \
--input DRB1-3123.fa.gz \
--n_haplotypes 12 \
--outdir testout \
-with-singularity  $PANGENOME_PATH/pangenome.img



```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=12 --mem=72g --gres=lscratch:10 pangenome.sh
```







