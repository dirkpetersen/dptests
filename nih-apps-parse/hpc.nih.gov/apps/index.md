

document.querySelector('title').textContent = 'Scientific Applications on NIH HPC Systems';Scientific Applications on NIH HPC Systems


|  |  |
| --- | --- |
| 
Application Areas
[Computational Chemistry](#compchem)
[Deep Learning](#deeplearning)
[Development Tools](/development/)
[High-Throughput Sequencing](#hts)
[Image Analysis](#image)
[Linkage/Phylogenetics](#linkage)
[Mass Spectrometry](#massspec)
[Mathematical/Statistics](#math)
[Molecular Modeling/Graphics](#model)
[Sequence Analysis](#seqanal)
[Structural Biology](#structbio)
[Systems Biology](#sysbiol)
[Utilities](#utils)
[Workflow Managers](#workflow)
 | 
The NIH HPC staff maintains several hundred scientific programs, packages and databases for our users.
Below is a list of system-installed software available on Biowulf and Helix. Click on the application name to
get to site-specific instructions on how to run a given package on the cluster, including links to the 
original application documentation.
 

In almost all cases, applications are made available through the use of
**[environment modules](modules.html)**.
 

Users are welcome to install additional applications in their own /home or /data areas, and create their own [personal modules](https://hpc.nih.gov/apps/modules.html#personal) for those applications. Some applications are difficult to install or require admin privileges; in those cases, you can email staff@hpc.nih.gov with a request to install the application.
  |


Computational Chemistry
[back to top](/apps/)  


[Acemd](acemd.html) (3.5.1) ACEMD is a high performance molecular dynamics code for biomolecular systems designed specifically for NVIDIA GPUs. Simple and fast, ACEMD uses very similar commands and input files of NAMD and output files as NAMD or Gromacs. 


[AMBER](AMBER.html) (22) AMBER (Assisted Model Building with Energy Refinement) is a package of molecular simulation programs.


[AMPL](ampl.html) (1.3.0) The Accelerating Therapeutics for Opportunites in Medicine (ATOM) Consortium Modeling PipeLine for Drug Discovery. AMPL is an open-source, modular, extensible software pipeline for building and sharing models to advance in silico drug discovery.


[APBS](APBS.html) (3.4.1) APBS (Adaptive Poisson-Boltzmann Solver) is a software package for the numerical solution of the Poisson-Boltzmann equation (PBE), one of the most popular continuum models for describing electrostatic interactions between molecular solutes in salty, aqueous media.


Autodock (4.2.6) Autodock is a suite of automated docking tools. It is designed to predict how small molecules, such as substrates or drug candidates, bind to a receptor of known 3D structure.


[AutodockVina](AutodockVina.html) (1\_1\_2) AutoDock Vina is a program for drug discovery, molecular docking and virtual screening, offering multi-core capability, high performance and enhanced accuracy and ease of use. It is closely tied to Autodock.


[CHARMM](charmm/index.html) (c46b1) CHARMM is a general and flexible software application for modeling the structure and behavior of molecular systems.


[cryoDRGN](cryoDRGN.html) (1.1.0) CryoDRGN is an algorithm that leverages the representation power of deep neural networks
to directly reconstruct continuous distributions of 3D density maps
and map per-particle heterogeneity of single-particle cryo-EM datasets.
It contains interactive tools to visualize a dataset’s distribution of per-particle
variability, generate density maps for exploratory analysis, extract particle
subsets for use with other tools and generate trajectories to visualize molecular motions.


[GAMESS](GAMESS.html) (30Sep22-R2) GAMESS is a general ab initio quantum chemistry package.


[Gaussian](Gaussian.html) (G16-C02) Gaussian is a connected system of programs for performing semiempirical and ab initio molecular orbital (MO) calculations.


[gromacs](Gromacs.html) (2022.4) Gromacs is a versatile package to perform molecular dynamics, i.e. simulate the Newtonian equations of motion for systems with hundreds to millions of particles. It is primarily designed for biochemical molecules like proteins and lipids that have a lot of complicated bonded interactions, but since GROMACS is extremely fast at calculating the nonbonded interactions (that usually dominate simulations) many groups are also using it for research on non-biological systems, e.g. polymers.



[NAMD](NAMD.html) (2.14) NAMD is a parallel molecular dynamics program for UNIX platforms designed for high-performance simulations in structural biology. VMD, the associated molecular visualization program, is also available.


[parmed](parmed.html) (4.0.0) Cross-program parameter and topology file editor and molecular mechanical simulator engine.


[Psi4](https://hpc.nih.gov/apps/psi4.html) (1.6.1) Psi4 is an ab-initio electronic structure code that supports various methods for calculating energies and gradients of molecular systems.


[Q-Chem](qchem.html) (5.0.1) Q-Chem is a comprehensive ab initio quantum chemistry package for accurate predictions of molecular structures, reactivities, and vibrational, electronic and NMR spectra.


[Schrodinger](schrodinger/index.html) (2023.1) A limited number of Schrödinger applications are available on the Biowulf cluster through the Molecular Modeling Interest Group. Most are available through the Maestro GUI.


[vmd](/apps/vmd.html) (1.9.3) VMD is a molecular visualization program for displaying, animating, and analyzing large biomolecular systems using 3-D graphics and built-in scripting. To use, type vmd at the prompt.


Deep Learning
[back to top](/apps/)  


[B-SOID](b-soid.html) (1.3) B-SOID (Behavioral Segmentation in Deeplabcut) is an 
unsupervised learning algorithm that serves to discover and classify behaviors that are not pre-defined by users.
It segregates statistically different, sub-second rodent behaviors with a single bottom-up perspective video cameraR by performing a novel expectation maximization
fitting of Gaussian mixture models on t-Distributed Stochastic Neighbor Embedding (t-SNE).


[BioBERT](BioBERT.html) (v20200409) BioBERT is a biomedical language representation model
designed for biomedical text mining tasks
such as biomedical named entity recognition, relation extraction, question answering, etc.


[BioGANs](biogans.html) (20210916) BioGANs is a novel application of Generative Adversarial Networks (GAN) to the synthesis of cells
imaged by fluorescence microscopy. It allows to infer the correlation between the spatial pattern
of different fluorescent proteins that reflects important biological functions.
The synthesized images capture these relationships, which are relevant for biological applications.


[cryoDRGN](cryoDRGN.html) (1.1.0) CryoDRGN is an algorithm that leverages the representation power of deep neural networks
to directly reconstruct continuous distributions of 3D density maps
and map per-particle heterogeneity of single-particle cryo-EM datasets.
It contains interactive tools to visualize a dataset’s distribution of per-particle
variability, generate density maps for exploratory analysis, extract particle
subsets for use with other tools and generate trajectories to visualize molecular motions.


[cuDNN](https://hpc.nih.gov/development/cuDNN.html) (8.0.3) The NVIDIA CUDA Deep Neural Network library (cuDNN) is a GPU-accelerated library of primitives for deep neural networks. cuDNN provides highly tuned implementations for standard routines such as forward and backward convolution, pooling, normalization, and activation layers.


[DanQ](DanQ.html) (20220825) DanQ is a hybrid convolutional and recurrent deep
neural network for quantifying the function of DNA sequences


[deepcadrt](deepcadrt.html) (0.1.0) Real-time denoising of fluorescence time-lapse imaging using deep self-supervised learning


[DeepCell-tf](deepcell-tf.html) (0.12.6) The DeepCell-tf library allows users to apply pre-existing models to imaging data
as well as to develop new deep learning models for single-cell analysis.
The library specializes in models for cell segmentation (whole-cell and nuclear) in 2D and 3D images
as well as cell tracking in 2D time-lapse datasets.
The models are applicable to data ranging from multiplexed images of tissues to dynamic live-cell imaging movies.


[DeepEmhancer](DeepEMhancer.html) (0.15) A deep learning solution for cryo-EM map.


[DeepLabCut](https://hpc.nih.gov/apps/DeepLabCut.html) (2.3.4) DeepLabCut is an open source toolbox
that builds on a state-of-the-art human pose estimation algorithm.
It allows training of a deep neural network by using limited training data
to precisely track user-defined features, so that the human labeling accuracy will be matched.


[deepmedic](deepmedic.html) (0.8.4) This project aims to offer easy access to Deep Learning for segmentation of structures of interest in biomedical 3D scans. It is a system that allows the easy creation of a 3D Convolutional Neural Network, which can be trained to detect and segment structures if corresponding ground truth labels are provided for training. The system processes NIFTI images, making its use straightforward for many biomedical tasks.


[DeepMM](DeepMM.html) (20220830) DeepMM implements fully automated de novo structure modeling method, MAINMAST, which builds three-dimensional models of a protein from a near-atomic resolution EM map. The method directly traces the protein’s main-chain and identifies Cα positions as tree-graph structures in the EM map.


[GCN\_Cancer](GCN_Cancer.html) (20221105) The GCN\_Cancer application employs graph convolutional network (GCN) models
to classify the gene expression data samples from The Cancer Genonme Atlas (TCAG)
as 33 designated tumor types or as normal. It has been trained on 10,340
cancer samples and 731 normal tissue samples from TCGA dataset.


[IsoNet](IsoNet.html) (0.2.1) IsoNet is a deep learning-based software package that iteratively reconstructs the missing-wedge information and increases signal-to-noise ratio, using the knowledge learned from raw tomograms. Without the need for sub-tomogram averaging, IsoNet generates tomograms with significantly reduced resolution anisotropy.


[keras](https://keras.io/) (2.4.3; 2.9.1) Keras is a deep learning API written in Python, running on top of the machine learning platform TensorFlow. It was developed with a focus on enabling fast experimentation. Being able to go from idea to result as fast as possible is key to doing good research.


[PyTom](PyTom.html) (1.0) PyTom is a software package for the analysis of volumetric data
obtained by cryo electron tomography (cryo-ET).
It covers a complete pipeline of processing steps
for tomogram reconstruction, localization of macromolecular complexes in tomograms,
fine alignment of subtomograms extracted at these locations, and their classification.


[ReLeaSE](ReLeaSE.html) (20220825) ReLeaSE (Reinforcement Learning for Structural Evolution) is an application for de-novo Drug Design
based on Reinforcement Learning. It integrates two deep neural networks: generative and
predictive, that are trained separately but are used jointly to generate novel targeted chemical libraries.
ReLeaSE uses simple representation of molecules by their simplified molecular input line entry specification (SMILES) strings only.


[RFdiffusion](RFdiffusion.html) (1.1.0) Rosetta Fold (RF) dissusion is an open source method for structure generation, with or without conditional information
(a motif, target etc). It can perform motif scaffolding, unconditional protein generation, and other tasks.


[SpliceAI](SpliceAI.html) (1.3.1) SpliceAI is a deep neural network that accurately predicts
splice junctions from an arbitrary pre-mRNA
transcript sequence, enabling precise prediction of
noncoding genetic variants that cause cryptic
splicing.



[talos](talos.html) (1.0) Talos is a hyperparameter optimization package for deep learning. It works with any Keras, TensorFlow (tf.keras) or PyTorch model, takes minutes to implement,
involves no new syntax to learn and adds zero new overhead to your workflow.


[TensorQTL](TensorQTL.html) (1.0.7) ensoorQTL leverages general-purpose libraries and graphics processing units
(GPUs) to achieve high efficiency of computations at low costR. Using PyTorch or TensorFlow
it allows > 200-fold decreases in runtime and ~ 5–10-fold reductions in cost when
running on GPUs relative to CPUs.


[Tybalt](Tybalt.html) (202208265) Tybalt implements a Variational EutoEncoder (VAE), a deep neural network approach capable of generating meaningful latent spaces for image and text data. Tybalt has been trained on The Cancer Genome Atlas (TCGA) pan-cancer RNA-seq data and used to identify specific patterns in the VAE encoded features.



[UNet](UNet.html) (20220825) U-Net is an image segmentation tool. It relies on the strong use of data augmentation to use the available
annotated samples more efficiently. The architecture consists of a contracting
path to capture context and a symmetric expanding path that enables precise localization.


High-Throughput Sequencing
[back to top](/apps/)  


[abyss](abyss.html) (2.3.5) Abyss represents Assembly By Short Sequences - a de novo, parallel, paired-end sequence assembler. The parallel version is implemented using MPI and is capable of assembling larger genomes.


[AfterQC](AfterQC.html) (0.9.7) Automatic Filtering, Trimming, Error Removing and Quality Control for fastq data. 



[angsd](angsd.html) (0.940)  ANGSD is a software for analyzing next generation sequencing data. The software can handle a number of different input types from mapped reads to imputed genotype probabilities. Most methods take genotype uncertainty into account instead of basing the analysis on called genotypes. This is especially useful for low and medium depth data. The software is written in C++ and has been used on large sample sizes. 


[AnnotSV](AnnotSV.html) (3.3.1) AnnotSV is a program designed for annotating Structural Variations (SV). This tool compiles functionally, regulatory and clinically relevant information and aims at providing annotations useful to i) interpret SV potential pathogenicity and ii) filter out SV potential false positives.



[ANNOVAR](ANNOVAR.html) (2020-06-08) ANNOVAR is an efficient software tool to utilize update-to-date information to functionally annotate genetic variants detected from diverse genomes.


[ascatNgs](ascatNgs.html) (4.5.0) AscatNGS contains the Cancer Genome Projects workflow implementation of the ASCAT copy number algorithm for paired end sequencing.


[augustus](augustus.html) (3.4.0) AUGUSTUS is a program that predicts genes in eukaryotic genomic sequences. 


[autoreject](autoreject.html) (0.3.1) This is a library to automatically reject bad trials and repair bad sensors in magneto-/electroencephalography (M/EEG) data.


[bam2fastq](bam2fastq.html) (1.1.0) This tool is used to extract raw sequences (with qualities) from bam files. 


[bamliquidator](bamliquidator.html) (1.5.2) bamliquidator is a set of tools for analyzing the density of short DNA sequence read alignments in the BAM file format.


[bamreadcount](bamreadcount.html) (cram-v0.0.1) Bam-readcount generates metrics at single nucleotide positions.
There are number of metrics generated which can be useful for filtering out false positive calls.


[bamtools](bamtools.html) (2.5.2) BamTools provides a fast, flexible C++ API & toolkit for reading, writing, and manipulating BAM files. 


[bamUtil](bamUtil.html) (1.0.15) bamUtil is a repository that contains several programs that perform operations on SAM/BAM files. All of these programs are built into a single executable, bam.


[Bartender](bartender.html) (1.1) Bartender is an accurate clustering algorithm to detect barcodes
and their abundances from raw next-generation sequencing data.
In contrast with existing methods that cluster based on sequence similarity alone,
Bartender uses a modified two-sample proportion test that also considers cluster size. This modification
results in higher accuracy and lower rates of under- and over-clustering artifacts.


[basespace\_cli](basespace_cli.html) (1.5.2) Command line interface for Illumina's BaseSpace


[bazam](bazam.html) (1.0.1) A tool to extract paired reads in FASTQ format from coordinate sorted BAM files


[bbtools](bbtools.html) (38.96) An extensive set of bioinformatics tools including bbmap (short read aligner), bbnorm (kmer based normalization), dedupe (deduplication and clustering of unaligned reads), reformat (formatting and trimming reads) and many more.


[bcl-convert](bcl-convert.html) (4.1.5) The Illumina BCL Convert is a standalone local software app that converts the Binary Base Call (BCL) files produced by Illumina sequencing systems to FASTQ files. BCL Convert also provides adapter handling (through masking and trimming) and UMI trimming and produces metric outputs.


[bcl2fastq](bcl2fastq.html) (2.20) a tool to handle bcl conversion and demultiplexing


[bedops](bedops.html) (2.4.41) Bedops is a suite of tools to address common questions raised in genomic studies - mostly with regard to overlap and proximity relationships between data sets - BEDOPS aims to be scalable and flexible, facilitating the efficient and accurate analysis and management of large-scale genomic data.


[bedtools](bedtools.html) (2.30.0) The BEDTools utilities allow one to address common genomics tasks such finding feature overlaps and computing coverage. In addition, one can develop sophisticated pipelines that answer complicated research questions by "streaming" several BEDTools together.


[BETA](beta.html) (1.0.7) Binding and expression target analysis (BETA) is a software package that integrates ChIP-seq of TFs or chromatin regulators with differential gene expression data to infer direct target genes. The combination of ChIP-seq and transcriptome analysis is a compelling approach to unravel the regulation of gene expression.



[bigSCale2](https://hpc.nih.gov/apps/bigSCale2.html) (20191119) bigSCale is a complete framework for the analysis and visualization of single cell data. It allows to cluster, phenotype, perform pseudotime analysis, infer gene regulatory networks and reduce large datasets in smaller datasets with higher quality.



[bioawk](bioawk.html) (1.0) Regular awk with support for several common biological data formats, including optionally gzip'ed BED, GFF, SAM, VCF, FASTA/Q and TAB-delimited formats with column names.


[biobakery\_workflows](biobakery_workflows.html) (3.1) bioBakery is a meta’omic analysis environment and collection of individual software
tools with the capacity to process raw shotgun sequencing data into actionable microbial community
feature profiles, summary reports, and publication-ready figures.
It includes a collection of preconfigured analysis modules also joined into workflows for reproducibility.
Each individual module has been developed to perform a particular task, e.g. quantitative taxonomic
profiling or statistical analysis.


[biobambam2](https://hpc.nih.gov/apps/biobambam2.html) (2.0.185-release-20221211202123) Tools for early stage alignment file processing.


[biom-format](biom-format.html) (2.1.12) tool (and library) to manipulate Biological Observation Matrix (BIOM) Format files


[bismark](bismark.html) (0.23.1) Bismark is a program to map bisulfite treated sequencing reads to a genome of interest and perform methylation calls in a single step. The output can be easily imported into a genome viewer, such as SeqMonk, and enables a researcher to analyse the methylation levels of their samples straight away.


[bonito](bonito.html) (0.7.2) A PyTorch Basecaller for Oxford Nanopore Reads


[bowtie](bowtie.html) (1.3.1) bowtie is an ultrafast, memory-efficient short read aligner geared toward quickly aligning large sets of short DNA sequences (reads) to large genomes.


[bowtie2](bowtie2.html) (2.5.1) A version of bowtie that's particularly good at aligning reads of about 50 up to 100s or 1,000s of characters, and particularly good at aligning to relatively long (e.g. mammalian) genomes


[BRASS](BRASS.html) (6.3.4) BRASS analyses one or more related BAM files of paired-end sequencing to determine potential rearrangement breakpoints.



[breakdancer](breakdancer.html) (1.4.5) provides genome-wide detection of structural variants from next generation paired-end sequencing reads.


[breseq](breseq.html) (0.37.1) breseq is a computational pipeline for finding mutations relative to a reference sequence in short-read DNA re-sequencing data. It is intended for haploid microbial genomes (<20 Mb).


[bsmap](bsmap.html) (2.90) BSMAP is a short reads mapping software for bisulfite sequencing reads. Bisulfite treatment converts unmethylated Cytosines into Uracils (sequenced as Thymine) and leave methylated Cytosines unchanged, hence provides a way to study DNA cytosine methylation at single nucleotide resolution. BSMAP aligns the Ts in the reads to both Cs and Ts in the reference.



[busco](busco.html) (5.4.7) BUSCO completeness assessments employ sets of Benchmarking Universal Single-Copy Orthologs from OrthoDB (www.orthodb.org) to provide quantitative measures of the completeness of genome assemblies, annotated gene sets, and transcriptomes in terms of expected gene content. 


[bwa](bwa.html) (0.7.17) BWA is a fast light-weighted tool that aligns short sequences to a sequence database, such as the human reference genome.


[bwa-mem2](bwa-mem2.html) (2.2.1) The next version of the bwa-mem algorithm in bwa. 


[CADD](CADD.html) (1.6.post1) CADD (Combined Annotation Dependent Depletion) is a tool for scoring the deleteriousness of single nucleotide variants as well as insertion/deletions variants in the human genome. Currently, it supports the builds: GRCh37/hg19 and GRCh38/hg38.


[Canu](canu.html) (2.1) Canu is a fork of the Celera Assembler designed for high-noise single-molecule sequencing (such as the PacBio RSII or Oxford Nanopore MinION). Canu will correct the reads, then trim suspicious regions (such as remaining SMRTbell adapter), then assemble the corrected and cleaned reads into unitigs.




[Canvas](Canvas.html) (1.40) Canvas is a tool for calling copy number variants (CNVs) from human DNA sequencing data.


[ccbrpipeliner](https://github.com/CCBR) (4.0.2) CCBR Pipeliner provides access to a set of best-practices NGS pipelines developed, tested, and benchmarked by experts at CCBR and NCBR for Biowulf. Contact CCBR\_Pipeliner@mail.nih.gov with questions


[ceas](ceas.html) (1.0.2) Cis-regulatory Element Annotation System is a tool designed to characterize genome-wide protein-DNA interaction patterns from ChIP-chip and ChIP-Seq of both sharp and broad binding factors. 


[cellprofiler](cellprofiler.html) (4.2.5) An open-source application for biological image analysis


[cellranger](cellranger.html) (7.2.0) Cell Ranger is a set of analysis pipelines that processes Chromium single cell 3’ RNA-seq output to align reads, generate gene-cell matrices and perform clustering and gene expression analysis.


[cellranger-arc](cellranger-arc.html) (2.0.2) Cell Ranger ARC is a set of analysis pipelines that process Chromium Single Cell Multiome ATAC + Gene Expression sequencing data to generate a variety of analyses pertaining to gene expression, chromatin accessibility and their linkage. Furthermore, since the ATAC and gene expression measurements are on the very same cell, we are able to perform analyses that link chromatin accessibility and gene expression.


[cellranger-atac](cellranger-atac.html) (2.1.0) Cell Ranger ATAC is a set of analysis pipelines that process Chromium Single Cell ATAC data. 


[cellsnp-lite](cellsnp-lite.html) (1.2.2) Efficient genotyping bi-allelic SNPs on single cells


[cgpBattenberg](cgpBattenberg.html) (3.5.3) Detect subclonality and copy number in matched NGS data


[checkm2](checkm2.html) (1.0.2) Rapid assessment of genome bin quality using machine learning.


[CHESS](chess.html) (0.3.7) The CHESS (Comparison of Hi-C Experiments using Structural Similarity) application
implements an algorithm for the comparison of chromatin contact maps and automatic
differential feature extraction.


[ChromHMM](ChromHMM.html) (1.23) ChromHMM is software for learning and characterizing chromatin states. 


[circexplorer2](circexplorer2.html) (2.3.8) A combined strategy to identify circular RNAs (circRNAs and ciRNAs) 


[circos](circos.html) (0.69-9) Circos is a program for the generation of publication-quality, circularly composited renditions of genomic data and related annotations. Circos is particularly suited for visualizing alignments, conservation and intra and inter-chromosomal relationships. Also, Circos is useful to visualize any type of information that benefits from a circular layout. Thus, although it has been designed for the field of genomics, it is sufficiently flexible to be used in other data domains.


[Clair3](Clair3.html) (0.1-r6) Clair3 is a small variant caller for Illumina, PacBio and ONT long reads. Compare to PEPPER (r0.4), Clair3 (v0.1) shows a better SNP F1-score with ≤30-fold of ONT data (precisionFDA Truth Challenge V2), and a better Indel F1-score, while runs generally four times faster. 


[clark](clark.html) (1.2.6.1) A method based on a supervised sequence classification using discriminative k-mers


[cnvkit](cnvkit.html) (0.9.9) Copy number variant detection from targeted DNA sequencing


[cnvnator](cnvnator.html) (0.4.1) CNVnator is a tool for CNV discovery and genotyping from depth of read mapping. 


[cogentap](cogentap.html) (2.0.1) Cogent NGS Analysis Pipeline (CogentAP) is bioinformatic software for analyzing RNA-seq NGS data generated using various takara kits.


[combp](combp.html) (0.50.6) A library to combine, analyze, group and correct p-values in BED files. Unique tools involve correction for spatial autocorrelation. This is useful for ChIP-Seq probes and Tiling arrays, or any data with spatial correlation.


[conpair](https://hpc.nih.gov/apps/conpair.html) (0.2) Concordance and contamination estimator for tumor–normal pairs



[crispresso](crispresso.html) (2.2.8) Software pipeline for the analysis of CRISPR-Cas9 genome editing outcomes from sequencing data



[crossmap](crossmap.html) (0.6.5) CrossMap is a program for convenient conversion of genome coordinates between different assemblies (e.g. mm9->mm10). It can convert SAM, BAM, bed, GTF, GFF, wig/bigWig, and VCF files


[csvkit](csvkit.html) (1.3.0) csvkit is a suite of command-line tools for converting to and working with CSV, the king of tabular file formats.


[cufflinks](cufflinks.html) (2.2.1) Cufflinks assembles transcripts and estimates their abundances in RNA-Seq samples. It accepts aligned RNA-Seq reads and assembles the alignments into a parsimonious set of transcripts. Cufflinks then estimates the relative abundances of these transcripts based on how many reads support each one.



[cutadapt](cutadapt.html) (4.4) cutadapt removes adapter sequences from DNA high-throughput
sequencing data. This is usually necessary when the read length of the
machine is longer than the molecule that is sequenced, such as in
microRNA data.


[cutruntools2](cutruntools2.html) (2.0) cutruntools2 is a major update of CutRunTools, including a set of new features specially designed for CUT&RUN and CUT&Tag experiments. Both of the bulk and single-cell data can be processed, analyzed and interpreted.


[DANPOS](DANPOS.html) (3.1.1) A toolkit for Dynamic Analysis of Nucleosome and Protein Occupancy by Sequencing, version 2


[DeepCell-tf](deepcell-tf.html) (0.12.6) The DeepCell-tf library allows users to apply pre-existing models to imaging data
as well as to develop new deep learning models for single-cell analysis.
The library specializes in models for cell segmentation (whole-cell and nuclear) in 2D and 3D images
as well as cell tracking in 2D time-lapse datasets.
The models are applicable to data ranging from multiplexed images of tissues to dynamic live-cell imaging movies.


[deepconsensus](deepconsensus.html) (0.3.1) DeepConsensus uses gap-aware sequence transformers to correct errors in Pacific Biosciences (PacBio) Circular Consensus Sequencing (CCS) data


[deepsignal](deepsignal.html) (2-0.1.3) A deep-learning method for detecting DNA methylation state from Oxford Nanopore sequencing reads.


[deeptools](deeptools.html) (3.5.4) deepTools is a suite of user-friendly tools for the visualization, quality control and normalization of data from deep-sequencing DNA sequencing experiments.



[deepvariant](deepvariant.html) (1.5.0) DeepVariant is an analysis pipeline that uses a deep neural network to call genetic variants from next-generation DNA sequencing data. 


[defuse](defuse.html) (0.8.1) deFuse is a software package for gene fusion discovery using RNA-Seq data. The software uses clusters of discordant paired end alignments to inform a split read alignment analysis for finding fusion boundaries. The software also employs a number of heuristic filters in an attempt to reduce the number of false positives and produces a fully annotated output for each predicted fusion.


[delly](delly.html) (1.1.6) DELLY is an integrated structural variant prediction method that can detect deletions, tandem duplications, inversions and translocations at single-nucleotide resolution in short-read massively parallel sequencing data. It uses paired-ends and split-reads to sensitively and accurately delineate genomic rearrangements throughout the genome.



[DETONATE](DETONATE.html) (1.11) DETONATE is a tool for evaluation of de novo transcriptome
assemblies from RNA-Seq data. It consists of two component packages, RSEM-EVAL and REF-EVAL. RSEM-EVAL is a reference-free evaluation method based on a novel probabilistic model that depends only on an assembly and the RNA-Seq reads used for its construction. REF-EVAL is a toolkit of reference-based measures.


[distiller-nf](distiller-nf.html) (0.3.3) A modular Hi-C mapping pipeline for reproducible data analysis, it was used for Micro-C analysis too. 


[dorado](dorado.html) (0.4.2) Dorado is a high-performance, easy-to-use, open source basecaller for Oxford Nanopore reads.


[drop](drop.html) (1.3.3) A pipeline to find aberrant events in RNA-Seq data.


[ea-utils](https://expressionanalysis.github.io/ea-utils/) (1.04.807) Command-line tools for processing biological sequencing data. Barcode demultiplexing, adapter trimming, etc. 


[edd](edd.html) (1.1.19) EDD is a ChIP-seq peak caller for detection of megabase domains of enrichment.


[encode-atac-seq-pipeline](encode-atac-seq-pipeline.html) (2.1.0) This pipeline is designed for automated end-to-end quality control and processing of ATAC-seq or DNase-seq data.


[EPACTS](EPACTS.html) (3.4.2) EPACTS (Efficient and Parallelizable Association Container Toolbox) is a versatile software pipeline to perform various statistical tests for identifying genome-wide association from sequence data through a user-friendly interface, both to scientific analysts and to method developers.



[eukdetect](https://hpc.nih.gov/apps/eukdetect.html) (1.2) EukDetect: Detect eukaryotes from shotgun metagenomic data


[exomiser](exomiser.html) (13.0.1) The Exomiser is a Java program that functionally annotates variants from whole-exome sequencing data starting from a VCF file.


[express](express.html) (1.5.1) eXpress is a streaming tool for quantifying the abundances of a set of target sequences from sampled subsequences. 


[fanc](fanc.html) (0.9.21) FAN-C is a toolkit for the analysis and visualization of Hi-C data. Beyond objects generated within FAN-C, the toolkit is largely compatible with Hi-C files from Cooler and Juicer.


[fastqc](fastqc.html) (0.11.9) It provide quality control functions to next gen sequencing data.


[fastqtools](fastqtools.html) (0.8.3) fastq-tools a collection of small and efficient programs for performing some common and uncommon tasks with FASTQ files. 


[fastq\_screen](fastq_screen.html) (0.15.3) FastQ Screen allows you to screen a library of sequences in FastQ format against a set of sequence databases so you can see if the composition of the library matches with what you expect.


[fastxtoolkit](fastxtoolkit.html) (0.0.14) The FASTX-Toolkit is a collection of command line tools for Short-Reads FASTA/FASTQ files preprocessing.


[fcs](https://hpc.nih.gov/apps/fcs.html) (0.4.0) FCS is a toolset to remove contaminant sequences from a genome assembly.


fgbio (2.0.2) The Fulcrum Genomics tools are a set of utilities for working with BAM files, VCF files, and Unique Molecular IDs. Theey are accessed as subprograms from a Java jar, like GATK or Picard.


fithic (2.0.8) Fit-Hi-C is a tool for assigning statistical confidence estimates to intra-chromosomal contact maps produced by genome-wide genome architecture assays such as Hi-C. 


[flanker](flanker.html) (0.1.5) Gene-flank analysis tool


[flashpca](flashpca.html) (2.0) FlashPCA performs fast principal component analysis (PCA) of single nucleotide polymorphism (SNP) data, similar to smartpca from EIGENSOFT (http://www.hsph.harvard.edu/alkes-price/software/) and shellfish (https://github.com/dandavison/shellfish). FlashPCA is based on the https://github.com/yixuan/spectra/ library.


[Flexbar](flexbar.html) (3.5.0) Flexbar preprocesses high-throughput sequencing data efficiently. It demultiplexes barcoded runs and removes adapter sequences. Moreover, trimming and filtering features are provided. Flexbar increases read mapping rates and improves genome and transcriptome assemblies. It supports next-generation sequencing data in fasta and fastq format, e.g. from Illumina and the Roche 454 platform


[flye](flye.html) (2.9.1) Fast and accurate de novo assembler for single molecule sequencing reads


[fqtools](fqtools.html) (2.0) Tools for manipulating fastq files


[freebayes](freebayes.html) (1.3.5) Bayesian haplotype-based polymorphism discovery and genotyping


[freec](freec.html) (11.6) Control-FREEC is a tool for detection of copy-number changes and allelic imbalances (including LOH) using deep-sequencing data 


[fuseq-wes](fuseq-wes.html) (1.0.0) Fusion Gene Detection Using Whole-Exome Sequencing Data in Cancer Patients 


[fusioncatcher](fusioncatcher.html) (1.33) FusionCatcher searches for novel/known somatic fusion genes, translocations, and chimeras in RNA-seq data (paired-end or single-end reads from Illumina NGS platforms like Solexa/HiSeq/NextSeq/MiSeq) from diseased samples.


[GATK](GATK.html) (4.3.0.0) GATK, from the Broad Institute, is a structured software library that makes writing efficient analysis tools using next-generation sequencing data very easy, and second it's a suite of tools for working with human medical resequencing projects such as 1000 Genomes and The Cancer Genome Atlas. These tools include things like a depth of coverage analyzers, a quality score recalibrator, a SNP/indel caller and a local realigner.


[gem](gem.html) (3.4) High resolution peak calling and motif discovery for ChIP-seq and ChIP-exo data


[Gemini](gemini.html) (0.30.2) GEMINI (GEnome MINIng) is designed to be a flexible framework for exploring genetic variation in the context of the wealth of genome annotations available for the human genome. By placing genetic variants, sample genotypes, and useful genome annotations into an integrated database framework, GEMINI provides a simple, flexible, yet very powerful system for exploring genetic variation for for disease and population genetics.


[GEMMA](GEMMA.html) (0.98.5) GEMMA is the software implementing the Genome-wide Efficient Mixed Model Association algorithm for a standard linear mixed model and some of its close relatives for genome-wide association studies (GWAS). 


[Genome Browser](Genome_Browser.html) (451) The [Genome Browser Mirror Fragments](https://hpcnihapps.cit.nih.gov/genome) is a mirror of the UCSC Genome Browser. The URL is https://hpcnihapps.cit.nih.gov/genome. Users can also access the MySQL databases, supporting files directly, and a huge number of associated executables.


genrich (0.6) Genrich is a peak-caller for genomic enrichment assays (e.g. ChIP-seq, ATAC-seq). It analyzes alignment files generated following the assay and produces a file detailing peaks of significant enrichment.


[GeoMX NGS Pipeline](geomx_ngs_pipeline.html) (2.3.3.10) The GeoMx NGS Pipeline, processes RNA-sequencing files (FASTQ files) from Illumina sequencers according to parameters defined in the Configuration File (which is generated from the GeoMx DSP run). The Pipeline processes information from these files and outputs .dcc files, which can then be uploaded to the GeoMx DSP system for data analysis.


[Gistic](gistic.html) (2.0.23) Facilitates sensitive and confident localization of the targets of focal somatic copy-number alteration in human cancers.


[gmap-gsnap](gmap-gsnap.html) (2021-12-17) A Genomic Mapping and Alignment Programs


[gopeaks](gopeaks.html) (1.0.0) GoPeaks is a peak caller designed for CUT&TAG/CUT&RUN sequencing data. GoPeaks by default works best with narrow peaks such as H3K4me3 and transcription factors. However, broad epigenetic marks like H3K27Ac/H3K4me1 require different the step, slide, and minwidth parameters.


[graphaligner](graphaligner.html) (1.0.17b) Seed-and-extend program for aligning long error-prone reads to genome graphs.


[gridss](https://hpc.nih.gov/apps/gridss.html) (2.13.2) GRIDSS is a module software suite containing tools useful for the detection of genomic rearrangements. GRIDSS includes a genome-wide break-end assembler, as well as a structural variation caller for Illumina sequencing data. GRIDSS calls variants based on alignment-guided positional de Bruijn graph genome-wide break-end assembly, split read, and read pair evidence.


[gtex\_rnaseq](gtex_rnaseq.html) (V10) This module makes available the tools used in the GTEX RNA-Seq pipeline.


gtool (0.7.5)  GTOOL is a program for transforming sets of genotype data for use with the programs SNPTEST and IMPUTE.

GTOOL can be used to (a) generate subsets of genotype data,
 (b) to convert genotype data between the PED file format and the FILE FORMAT used by SNPTEST and IMPUTE.


[gunc](https://hpc.nih.gov/apps/gunc.html) (1.0.5) Genome UNClutterer (GUNC) is a tool for detection of chimerism and contamination in prokaryotic genomes resulting from mis-binning of genomic contigs from unrelated lineages


[guppy](guppy.html) (6.5.7) Local accelerated basecalling for Nanopore data


[Hail](Hail.html) (0.2.99) Hail is an open-source, scalable framework for exploring and analyzing genomic data.


[hap.py](hap.py.html) (0.3.14) A set of programs based on htslib to benchmark variant calls against gold standard truth datasets.


[hicexplorer](https://hpc.nih.gov/apps/hicexplorer.html) (3.7.2) Tools to process, normalize and visualize Hi-C data


[hichipper](hichipper.html) (0.7.7) hichipper is a preprocessing and QC pipeline for HiChIP data. This package takes output from a HiC-Pro run and a sample manifest file (.yaml) that coordinates optional high-quality peaks (identified through ChIP-Seq) and restriction fragment locations (see folder here) as input and produces output that can be used to 1) determine library quality, 2) identify and characterize DNA loops and 3) interactively visualize loops.


[hicpro](hicpro.html) (3.1.0) HiC-Pro: An optimized and flexible pipeline for Hi-C data processing


hic\_breakfinder (1.0) A framework that integrates optical mapping, high-throughput chromosome conformation capture (Hi-C), and whole genome sequencing to systematically detect SVs in a variety of normal or cancer samples and cell lines.


[hifiasm](hifiasm.html) (0.19.5)  Hifiasm is a fast haplotype-resolved de novo assembler initially designed for PacBio HiFi reads.
 Its latest release supports telomere-to-telomere assembly by utilizing ultralong Oxford Nanopore reads.
 It can produce better haplotype-resolved assemblies when given parental short reads or Hi-C data.


[hint](hint.html) (2.2.7) a computational method to detect CNVs and Translocations from Hi-C data. 


[hipstr](hipstr.html) (0.7) Tool for genotyping short tandom repeats from Illumina sequencing data


[hisat](hisat.html) (2.2.2.1-ngs3.0.1) HISAT is a fast and sensitive spliced alignment program which uses Hierarchical Indexing for Spliced Alignment of Transcripts.


[HLA-LA](HLA-LA.html) (1.0.3) Fast HLA type inference from whole-genome data. Previously known as HLA-PRG-LA.


[HMMRATAC](HMMRATAC.html) (1.2.10) HMMRATAC peak caller for ATAC-seq data


[homer](homer.html) (4.11.1) HOMER (Hypergeometric Optimization of Motif EnRichment) is a suite of tools for Motif Discovery and ChIP-Seq analysis.


[htgts](htgts.html) (2) High-Throughput Genome-Wide Translocation Sequencing pipeline


[htsbox](htsbox.html) (r346) HTSbox is a fork of early HTSlib. It is a collection of small experimental tools manipulating HTS-related files. While some of these tools are already part of the official SAMtools package, others are for niche use cases.


[htseq](htseq.html) (2.0.4) HTSeq is a Python package that provides infrastructure to process data from high-throughput sequencing assays.


[humann](humann.html) (3.6.0) HUMAnN is a pipeline for efficiently and accurately profiling the presence/absence and abundance of microbial pathways in a community from metagenomic or metatranscriptomic sequencing data (typically millions of short DNA/RNA reads).


[IDBA](idba.html) (1.1.3) IDBA is a practical iterative De Bruijn Graph De Novo Assembler for sequence assembly in bioinfomatics. Most assemblers based on de Bruijn graph build a de Bruijn graph with a specific k to perform the assembling task. For all of them, it is very crucial to find a specific value of k. If k is too large, there will be a lot of gap problems in the graph. If k is too small, there will a lot of branch problems. IDBA uses not only one specific k but a range of k values to build the iterative de Bruijn graph. It can keep all the information in graphs with different k values.


[iDiffIR](iDiffIR.html) (20220121) iDiffIR is a tool for identifying differential IR from RNA-seq data.
It accepts any sorted, indexed BAM file for single- or paired-end reads.


[IDR](idr.html) (2.0.3) The IDR (Irreproducible Discovery Rate) framework is a uniﬁed approach to measure the reproducibility of ﬁndings identiﬁed from replicate experiments and provide highly stable thresholds based on reproducibility. The IDR method compares a pair of ranked lists of identifications (such as ChIP-seq peaks).


[IGV](IGV.html) (2.12.3) The Integrative Genomics Viewer is a high-performance visualization tool for interactive exploration of large, integrated genomic datasets. 


[IGVTools](IGV.html) (2.12.3) IGVTools provides utilities for working with ascii file formats used by the Integrated Genome Viewer. The files can be sorted, tiled, indexed, and counted.


[IMPUTE](IMPUTE.html) (2.3.2) Impute is a program for estimating ("imputing") unobserved genotypes in SNP association studies.


[InterVar](InterVar.html) (2.1.2, 2.1.3) In 2015, the American College of Medical Genetics and Genomics (ACMG) and the Association for Molecular Pathology (AMP) published
updated standards and guidelines for the clinical interpretation of sequence variants with respect to human diseases on the basis
of 28 criteria. However, variability between individual interpreters can be extensive because of reasons such as the different understandings of these guidelines and the lack of standard algorithms for implementing them, yet computational tools for semi-automated variant interpretation are not available. To address these problems, InterVar implements these criteria to help human reviewers interpret the clinical significance of variants. InterVar can take a pre-annotated or VCF file as input and generate automated interpretation on 18 criteria.


intervene (0.6.5) a tool for intersection and visualization of multiple genomic region sets


[iva](iva.html) (1.0.11) IVA is a de novo assembler designed to assemble virus genomes that have no repeat sequences, using Illumina read pairs sequenced from mixed populations at extremely high and variable depth.


[iVar](iVar.html) (1.3.1) iVar is a computational package that contains functions broadly useful for viral amplicon-based sequencing. Additional tools for metagenomic sequencing are actively being incorporated into iVar. While each of these functions can be accomplished using existing tools, iVar contains an intersection of functionality from multiple tools that are required to call iSNVs and consensus sequences from viral sequencing data across multiple replicates.


[Juicer](juicer.html) (1.6) A One-Click System for Analyzing Loop-Resolution Hi-C Experiments


[jvarkit](jvarkit.html) (20211013) Java tools for bioinformatics


[KAT](kat.html) (2.4.2)  KAT (K-mer Analysis Toolkit) is a suite of tools that analyse Jellyfish hashes or sequence files (fasta or fastq) using kmer counts.


[kb-python](kb-python.html) (0.27.3) kb-python is a python package for processing single-cell RNA-sequencing. It wraps the kallisto | bustools single-cell RNA-seq command line tools in order to unify multiple processing workflows. 


[KMC](kmc.html) (3.1.0) KMC is a disk-based programm for counting k-mers from (possibly gzipped) FASTQ/FASTA files


[KmerGO2](https://hpc.nih.gov/apps/KmerGO2.html) (2.0.1) KmerGO is a user-friendly tool to identify the group-specific sequences on two groups of high throughput sequencing datasets. 


[kneaddata](kneaddata.html) (0.11.0) KneadData is a tool designed to perform quality control on metagenomic and metatranscriptomic sequencing data, especially data from microbiome experiments.


[kraken](kraken.html) (2.1.2) Kraken is a system for assigning taxonomic labels to short DNA sequences, usually obtained through metagenomic studies


[lefse](lefse.html) (1.0.8) LEfSe (Linear discriminant analysis Effect Size) determines the features (organisms, clades, operational taxonomic units, genes, or functions) most likely to explain differences between classes by coupling standard tests for statistical significance with additional tests encoding biological consistency and effect relevance.


[LJA](LJA.html) (0.1) The La Jolla Assembler (LJA) is a tool for assemblies of long and accurate reads. It reduces the error rate in
these reads by three orders of magnitude (making them nearly error-free) and constructs the de Bruijn
graph for large genomes and large k-mer sizes. Since the de Bruijn graph constructed for a fixed k-mer
size is typically either too tangled or too fragmented, LJA uses a new concept of a multiplex de Bruijn
graph with varying k-mer sizes.


[locuszoom](locuszoom.html) (1.3) LocusZoom is designed to facilitate viewing of local association results together with useful information about a locus, such as the location and orientation of the genes it includes, linkage disequilibrium coefficients and local estimates of recombination rates


[lofreq](lofreq.html) (2.1.5) LoFreq is a fast and sensitive variant-caller for inferring SNVs and indels from next-generation sequencing data.


[lordec](lordec.html) (141) LoRDEC processes data coming from high throughput sequencing machines of the second and third generations. These data are called sequencing reads, or simply reads for short. Technically speaking it processes short reads and long reads to correct errors in the long reads.


[lumpy](lumpy.html) (0.3.1) A probabilistic framework for structural variant discovery.


[mach2qtl](mach2qtl.html) (1.1.3) mach2qtl uses dosages/posterior probabilities inferred with MACH as predictors in a linear regression to test association with a quantitative trait


[macs](macs.html) (3) Model-based Analysis of ChIP-Seq (MACS) on short reads sequencers such as Genome Analyzer (Illumina / Solexa). MACS empirically models the length of the sequenced ChIP fragments, which tends to be shorter than sonication or library construction size estimates, and uses it to improve the spatial resolution of predicted binding sites. MACS also uses a dynamic Poisson distribution to effectively capture local biases in the genome sequence, allowing for more sensitive and robust prediction.


[mafft](mafft.html) (7.475) Multiple alignment program for amino acid or nucleotide sequences



[MAGeCK](MAGeCK.html) (0.5.9.2) MAGeCK is Model-based Analysis of Genome-wide CRISPR/Cas9 Knockout (MAGeCK) method for prioritizing
single-guide RNAs, genes and pathways in genome-scale CRISPR/Cas9 knockout screens. It demonstrates better performance compared with other methods, identifies both positively and negatively selected genes
simultaneously, and reports robust results across different experimental conditions.



[mageck-vispr](mageck-vispr.html) (0.5.6) MAGeCK-VISPR is a comprehensive quality control, analysis and visualization workflow for CRISPR/Cas9 screens.


[magicblast](magicblast.html) (1.7.0) Magic-BLAST is a tool for mapping large next-generation RNA or DNA sequencing runs against a whole genome or transcriptome. Each alignment optimizes a composite score, taking into account simultaneously the two reads of a pair, and in case of RNA-seq, locating the candidate introns and adding up the score of all exons. This is very different from other versions of BLAST, where each exon is scored as a separate hit and read-pairing is ignored.


[MAJIQ](majiq.html) (2.4) Modeling Alternative Junction Inclusion Quantification. MAJIQ and Voila are two software packages that together define, quantify, and visualize local splicing variations (LSV) from RNA-Seq data.


[manorm](manorm.html) (1.1.4) MAnorm is for quantitative comparison of ChIP-Seq data sets describing transcription factor binding sites and epigenetic modifications. The quantitative binding differences inferred by MAnorm showed strong correlation with both the changes in expression of target genes and the binding of cell type-specific regulators.


[manta](manta.html) (1.6.0-fork-jykr) Structural variant and indel caller for mapped sequencing data 


maps (1.1.0) a set of multiple scripts used to analyze PLAC-Seq and HiChIP data.


[mash](mash.html) (2.3) mash is a command line tool and library to provide fast genome and metagenome distance estimation using MinHash. Only command line tool is installed


[mbg](mbg.html) (1.0.15) Minimizer based sparse de Bruijn Graph constructor.


[medaka](medaka.html) (1.10.0) medaka is a tool to create a consensus sequence from nanopore sequencing data. This task is performed using neural networks applied from a pileup of individual sequencing reads against a draft assembly.


[medusa](medusa.html) (1.6) A draft genome scaffolder that uses multiple reference genomes in a graph-based approach.


[megadepth](megadepth.html) (1.2.0) 
[MEGAHIT](megahit.html) (1.2.9) MEGAHIT is a single node assembler for large and complex metagenomics NGS reads, such as soil. It makes use of succinct de Bruijn graph (SdBG) to achieve low memory assembly. MEGAHIT can optionally utilize a CUDA-enabled GPU to accelerate its SdBG construction.


[megalodon](megalodon.html) (2.5.0) Megalodon is a research command line tool to extract high accuracy modified base and sequence variant calls from raw nanopore reads by anchoring the information rich basecalling neural network output to a reference genome/transcriptome.


[MEGAN](MEGAN.html) (6.24.11) MEtaGenome ANalyzer that takes a file of reads and a Blast output from comparison against a reference genome, and automatically calculate a taxonomic classification of the reads and if desired, a functional classification.


[merfin](merfin.html) (1.1) Improved variant filtering and polishing via k-mer validation


[merqury](merqury.html) (1.3) Evaluate genome assemblies with k-mers and more


[Meryl](https://hpc.nih.gov/apps/Meryl.html) (0.0) Meryl: a genomic k-mer counter (and sequence utility) with nice features. It is built into the Celera Assembler
and is also available as a stand-alone application.
Meryl uses a sorting-based approach that sorts the k-mers in lexicographical order.


[metabat](metabat.html) (2.15) MetaBAT: A robust statistical framework for reconstructing genomes from metagenomic data


[metal](metal.html) (2020-05-05) The METAL software is designed to facilitate meta-analysis of large datasets (such as several whole genome scans) in a convenient, rapid and memory efficient manner. 


[metaphlan](metaphlan.html) (4.0.3) MetaPhlAn is a computational tool for profiling the composition of microbial communities (Bacteria, Archaea, Eukaryotes and Viruses) from metagenomic shotgun sequencing data (i.e. not 16S) with species-level. With the newly added StrainPhlAn module, it is now possible to perform accurate strain-level microbial profiling.


miniasm (0.3.r179) Ultrafast de novo assembly for long noisy reads (though having no consensus step)


[minimac](minimac.html) (4 (1.0.1)) minimac is a low memory, computationally efficient implementation of the MaCH algorithm for genotype imputation. It is designed to work on phased genotypes and can handle very large reference panels with hundreds or thousands of haplotypes. 'mini' refers to the low amount of computational resources it needs. 


[minimap2](minimap2.html) (2.26) Minimap2 is a fast sequence mapping and alignment program that can find overlaps between long noisy reads, or map long reads or their assemblies to a reference genome optionally with detailed alignment (i.e. CIGAR).


[mirdeep2](mirdeep.html) (0.1.3) miRDeep2 is a completely overhauled tool which discovers microRNA genes by analyzing sequenced RNAs. 


[misopy](misopy.html) (0.5.4) MISO (Mixture-of-Isoforms) is a probabilistic framework that quantitates the expression level of alternatively spliced genes from RNA-Seq data


[mitosuite](mitosuite.html) (1.0.9b) mitosuite is a graphical tool for human mitochondrial genome profiling in massively parallel sequencing


[mixcr](mixcr.html) (4.3.2) MiXCR is a universal software for fast and accurate analysis of T- and B- cell receptor repertoire sequencing data.


[mixer](mixer.html) (1.3) MiXeR is Causal Mixture Model for GWAS summary statistics. The version(1.3) contains a Python port of MiXeR, wrapping the C/C++ core. Also data preprocessing code sumstats.py is included too.



[modbamtools](modbamtools.html) (0.4.8) A set of tools to manipulate and visualize DNA/RNA base modification data that are stored in bam format.


[modkit](modkit.html) (0.1.11) A bioinformatics tool for working with modified bases from Oxford Nanopore. Specifically for converting modBAM to bedMethyl files using best practices, but also manipulating modBAM files and generating summary statistics.


[modphred](modphred.html) (1.0c) modPhred is a pipeline for detection, annotation and visualisation of DNA/RNA modifications.


[mosaicforecast](mosaicforecast.html) (0.0.1) a machine learning method that leverages read-based phasing and read-level features to accurately detect mosaic SNVs (SNPs, small indels) from NGS data.It builds on existing algorithms to achieve a multifold increase in specificity. 


[mosdepth](mosdepth.html) (0.3.3) Fast BAM/CRAM depth calculation for WGS, exome, or targeted sequencing. 


[mtoolbox](mtoolbox.html) (1.2.1) A bioinformatics pipeline aimed at the analysis of mitochondrial DNA (mtDNA) in high throughput sequencing studies.


[multiqc](multiqc.html) (1.14) aggregates results for various frequently used bioinformatics tools across multiple samples into a nice visual report


[MuSE](http://bioinformatics.mdanderson.org/main/MuSE) (2.0.1) MuSE is an approach to somatic variant calling based on the F81 Markov substitution model for molecular evolution, which models the evolution of the reference allele to the allelic composition of the matched tumor and normal tissue at each genomic locus.


[muTect](muTect.html) (1.1.7) MuTect is a method developed at the Broad Institute for the reliable and accurate identification of somatic point mutations in next generation sequencing data of cancer genomes.


[MutSig](MutSig.html) (1.41) MutSig analyzes lists of mutations discovered in DNA sequencing, to identify genes that were mutated more often than expected by chance given background mutation processes.


[nanopolish](nanopolish.html) (0.14.0) nanopolish is a software package for signal-level analysis of Oxford Nanopore sequencing data. Nanopolish can calculate an improved consensus sequence for a draft genome assembly, detect base modifications, call SNPs and indels with respect to a reference genome and more (see Nanopolish modules, below).




[neusomatic](neusomatic.html) (0.2.1) NeuSomatic is based on deep convolutional neural networks for accurate somatic mutation detection. With properly trained models, it can robustly perform across sequencing platforms, strategies, and conditions. NeuSomatic summarizes and augments sequence alignments in a novel way and incorporates multi-dimensional features to capture variant signals effectively. It is not only a universal but also accurate somatic mutation detection method.


[NGMLR](ngmlr.html) (0.2.7) NGMLR is a long-read mapper designed to align PacBio or Oxford Nanopore (standard and ultra-long) to a reference genome with a focus on reads that span structural variations.


[ngsplot](ngsplot/) (2.63) ngsplot is an easy-to-use global visualization tool for next-generation sequencing data.


[novocraft](novocraft.html) (4.03.08) Package includes aligner for single-ended and paired-end reads from the Illumina Genome Analyser. Novoalign finds global optimum alignments using full Needleman-Wunsch algorithm with affine gap penalties.


Octopus (0.7.4) Octopus is a mapping-based variant caller that implements several calling models within a unified haplotype-aware framework. Octopus takes inspiration from particle filtering by constructing a tree of haplotypes and dynamically pruning and extending the tree based on haplotype posterior probabilities in a sequential manner. This allows octopus to implicitly consider all possible haplotypes at a given loci in reasonable time.


[ont-fast5-api](ont-fast5-api.html) (4.1.0) Tools to manipulate HDF5 files of the Oxford Nanopore .fast5 file format


[pairtools](pairtools.html) (0.3.0; 1.0.2) Pairtools is a simple and fast command-line framework to process sequencing data
from a Hi-C experiment. Pairtools perform various operations on Hi-C pairs and occupy the middle position in a typical Hi-C data processing pipeline. Pairtools aim to be an all-in-one tool for processing Hi-C pairs.


[panaroo](panaroo.html) (1.3.3) An updated pipeline for pangenome investigation 


[parabricks](parabricks.html) (4.0.0) The Clara Parabricks toolkit is a set of GPU-accelerated genome analysis tools for secondary analysis of next generation sequencing data.


[PartekFlow](https://partekflow.cit.nih.gov) (10.0.23.0720) Web interface designed specifically for the analysis needs of next generation sequencing applications including RNA, small RNA, and DNA sequencing.


[pb-cpg-tools](pb-cpg-tools.html) (2.3.1) Tools for analyzing CpG/5mC data from PacBio HiFi reads aligned to a reference genome


[pbipa](pbipa.html) (1.8.0) Improved Phased Assembler (IPA) is the official PacBio software for HiFi genome assembly. IPA was designed to utilize the accuracy of PacBio HiFi reads to produce high-quality phased genome assemblies. 


[pbsuite](pbsuite.html) (15.8.24) The PBSuite contains two projects created for analysis of Pacific Biosciences long-read sequencing data: PBHoney and PBJelly. PBHoney is an implementation of two variant-identification approaches designed to exploit the high mappability of long reads (i.e., greater than 10,000 bp). PBHoney considers both intra-read discordance and soft-clipped tails of long reads to identify structural variants. PBJelly is a highly automated pipeline that aligns long sequencing reads (such as PacBio RS reads or long 454 reads in fasta format) to high-confidence draft assembles. PBJelly fills or reduces as many captured gaps as possible to produce upgraded draft genomes.


[pear](pear.html) (0.9.11) PEAR is an ultrafast, memory-efficient and highly accurate pair-end read merger. It is fully parallelized and can run with as low as just a few kilobytes of memory.


[peddy](peddy.html) (0.4.8) peddy is used to compare sex and familial relationships given in a PED file with those inferred from a VCF file


[PennCNV](PennCNV.html) ( 1.0.5) PennCNV is a free software tool for Copy Number Variation (CNV) detection
from SNP genotyping arrays. Currently it can handle signal intensity data
from Illumina and Affymetrix arrays. With appropriate preparation of file format,
it can also handle other types of SNP arrays and oligonucleotide arrays.


[PEPATAC](PEPATAC.html) (0.10.3) PEPATAC is a robust pipeline for Assay for Transposase-Accessible Chromatin using sequencing (ATAC-seq) built on a loosely coupled modular framework.
It may be easily applied to ATAC-seq projects of any size,
from one-off experiments to large-scale sequencing projects.
It is optimized on unique features of ATAC-seq data to be fast and accurate
and provides several unique analytical approaches.


[picard](picard.html) (3.1.0) Picard comprises Java-based command-line utilities that manipulate SAM files, and a Java API (SAM-JDK) for creating new programs that read and write SAM files. Both SAM text format and SAM binary (BAM) format are supported.


[picrust](picrust.html) (2.5.2) PICRUSt is a bioinformatics software package designed to predict metagenome functional content from marker gene (e.g., 16S rRNA) surveys and full genomes.


[pindel](pindel.html) (0.2.5) Pindel can detect breakpoints of large deletions, medium sized insertions, inversions, tandem duplications and other structural variants at single-based resolution from next-gen sequence data. It uses a pattern growth approach to identify the breakpoints of these variants from paired-end short reads.



[Platypus](platypus.html) (0.8.1) tool for variant-detection in high-throughput sequencing data. 


[plink](plink.html) (3.6-alpha) PLINK is whole genome association analysis toolset, designed to perform a range of basic, large-scale analyses in a computationally efficient manner.




[pod5](pod5.html) (0.5.0) Tool for manipulating the pod5 format of nanopore reads


[pomoxis](pomoxis.html) (0.3.10) Pomoxis comprises a set of basic bioinformatic tools tailored to nanopore sequencing. Notably tools are included for generating and analysing draft assemblies. Many of these tools are used by the research data analysis group at Oxford Nanopore Technologies.


[porechop](porechop.html) (0.2.4) Trim/demultiplex Oxford Nanopore reads


[poretools](poretools.html) (0.6.1a1) Poretools is a toolkit for manipulating and exploring nanopore sequencing data sets. Poretools operates on individual FAST5 files, directory of FAST5 files, and tar archives of FAST5 files. 


[preseq](preseq.html) (3.1.2) predicting library complexity and genome coverage in high-throughput sequencing


[PRINSEQ](prinseq.html) (0.20.4) PRINSEQ is a tool that generates summary statistics of sequence and quality data and that is used to filter, reformat and trim next-generation sequence data. It is particular designed for 454/Roche data, but can also be used for other types of sequence data.


[proseq](proseq.html) (2.0) proseq-2.0 is a pipeline for preprocesses and alignment
of run-on sequencing (PRO/GRO/ChRO-seq) data
from Single-Read or Paired-End Illumina Sequencing



[pvactools](pvactools.html) (4.0.1) pVACtools is a cancer immunotherapy suite consisting of pVACseq, pVACfuse, pVACvector


[pychopper](pychopper.html) (2.7.1) Pychopper v2 is a tool to identify, orient and trim full-length Nanopore cDNA reads. The tool is also able to rescue fused reads.


[pyclone](pyclone.html) (0.13.1) PyClone is statistical model and software tool designed to infer the prevalence of point mutations in heterogeneous cancer samples. 


[pycoQC](pycoQC.html) (2.5.2) pycoQC is a new tool to generate interactive quality control metrics and plots
from basecalled nanopore reads
or summary files generated by the basecallers Albacore, Guppy or MinKNOW.
pycoQC has several novel features, including:
1) python support for creation of dynamic D3.js visualizations and interactive data exploration in Jupyter Notebooks;
2) simple command line interface to generate customizable interactive HTML reports; and
3) multiprocessing FAST5 feature extraction program to generate a summary file directly
from FAST5 files.


[pyega3](pyega3.html) (5.0.2) A download client for the European Genome-phenome Archive (EGA). The EGA is designed to be a repository for all types of sequence and genotype experiments, including case-control, population, and family studies.


qctool (2.2) QCTOOL is a command-line utility program for basic quality control of gwas datasets. 


[qualimap](qualimap.html) (2.2.1) A platform-independednt application written in Java and R that provides both a GUI and a co
mmand-line interface to facilitate the quality control of alignment sequencing data.


[quast](quast.html) (5.2.0) QUAST stands for QUality ASsessment Tool. The tool evaluates genome assemblies by computing various metrics. The package includes the general QUAST tool for genome assemblies, MetaQUAST, the extension for metagenomic datasets, and Icarus, interactive visualizer for these tools.


racon (1.4.3) Ultrafast consensus module for raw de novo genome assembly of long uncorrected reads.


[ragtag](ragtag.html) (2.1.0) RagTag is a collection of software tools for scaffolding and improving modern genome assemblies.


[raremetal](raremetal.html) (4.15.1) RAREMETAL is a computationally efficient tool for meta-analysis of rare variants using sequencing or genotyping array data.


[Rcorrector](Rcorrector.html) (1.0.5) Rcorrector implements a k-mer based method to correct random sequencing errors
in Illumina RNA-seq reads. Rcorrector uses a De Bruijn graph to compactly represent
all trusted k-mers in the input reads. Unlike WGS read correctors,
which use a global threshold to determine trusted k-mers, Rcorrector computes a local
threshold at every position in a read.



[regtools](https://hpc.nih.gov/apps/regtools.html) (1.0.0) regtools:Tools that integrate DNA-seq and RNA-seq data to help interpret mutations in a regulatory and splicing context


[RepeatMasker](repeatmasker.html) (4.1.5) RepeatMasker is a program that screens DNA sequences for interspersed repeats and low complexity DNA sequences. The output of the program is a detailed annotation of the repeats that are present in the query sequence as well as a modified version of the query sequence in which all the annotated repeats have been masked (default: replaced by Ns). On average, almost 50% of a human genomic DNA sequence currently will be masked by the program. 


[repeatmodeler](repeatmodeler.html) (2.0.1) RepeatModeler is a de novo transposable element (TE) family identification and modeling package. RepeatModeler assists in automating the runs of the various algorithms given a genomic database, clustering redundant results, refining and classifying the families and producing a high quality library of TE families suitable for use with RepeatMasker and ultimately for submission to the Dfam database (http://dfam.org).


[REViewer](REViewer.html) (0.2.7) REViewer is a tool for visualizing alignments of reads in regions containing tandem repeats. REViewer requires a BAMlet with graph-realigned reads generated by ExpansionHunter and the corresponding variant catalog.


[rgt](rgt.html) (0.13.2) Regulatory Genomics Toolbox: Python library and set of tools for the integrative analysis of high throughput regulatory genomics data. http://www.regulatory-genomics.org


[rilseq](rilseq.html) (0.82) RILseq computational protocol


[rnaseqc](rnaseqc.html) (2.4.2) RNA-SeQC is a java program which computes a series of quality control metrics for RNA-seq data. 


[rockhopper](rockhopper.html) (2.0.3) Rockhopper is a comprehensive and user-friendly system for computational analysis of bacterial RNA-seq data. As input, Rockhopper takes RNA sequencing reads output by high-throughput sequencing technology (FASTQ, QSEQ, FASTA, SAM, or BAM files)


[ROSE](https://hpc.nih.gov/apps/ROSE.html) (1.3.1) ROSE (Rank Ordering of Super-Enhancers) is tool for 
(1) creating stitched enhancers, and 
(2) separating super-enhancers from typical enhancers. 
given sequencing data (.bam) and a file of previously identified constituent enhancers (.gff)


[RSD](RSD.html) (1.1.7) Reciprocal Smallest Distance (RSD) is a pairwise orthology algorithm that uses global sequence alignment and maximum likelihood evolutionary distance between sequences to accurately detects orthologs between genomes.


[rsem](rsem.html) (1.3.3) RSEM is a software package for estimating gene and isoform expression levels from RNA-Seq data. 


[rseqc](rseqc.html) (5.0.3) Rseqc comprehensively evaluate RNA-seq datasets generated from clinical tissues or other well annotated organisms such as mouse, fly and yeast. 


[rtg-tools](rtg-tools.html) (3.8.4) variant detection for singletons, families, large pedigrees and populations, cancer, structural variant and CNV analysis, and microbial and metagenomic analysis


[rvtests](rvtests.html) (2.1.0) Rare Variant tests is a flexible software package for genetic association studies. It is designed to support unrelated individual or related (family-based) individuals


[sailfish](sailfish.html) (0.10.0) Sailfish is a tool for transcript quantification from RNA-seq data. It requires a set of target transcripts (either from a reference or de-novo assembly) to quantify. All that is needed to run sailfish is a fasta file containing your reference transcripts and a (set of) fasta/fastq file(s) containing your RNA-Seq reads.


[salmon](salmon.html) (1.10.0) a tool for quantifying the expression of transcripts using RNA-seq data. 


[SalmonTE](salmonte.html) (0.4) SalmonTE is an ultra-Fast and Scalable Quantification Pipeline of Transpose Element (TE) Abundances.


[sambamba](sambamba.html) (0.8.2) Sambamba is a high performance modern robust and fast tool (and library), written in the D programming language, for working with SAM and BAM files. Current parallelised functionality is an important subset of samtools functionality, including view, index, sort, markdup, and depth.


[samblaster](samblaster.html) (0.1.26) samblaster is a program for marking duplicates and finding discordant/split read pairs in read-id grouped paired-end SAM files. When marking duplicates, samblaster will use about 20MB per 1M read pairs. In a read-id grouped SAM file all alignments for a read-id (QNAME) are continuous. Aligners naturally produce such files. They can also be created by sorting a SAM file by read-id. 


[samtools](samtools.html) (1.17) The samtools package now provides samtools, bcftools, tabix, and the underlying htslib library. 


[scallop](scallop.html) (0.10.5) Scallop is a reference-based transcript assembler.


[scalpel](scalpel.html) (0.5.4) Bioinformatics pipeline for discovery of genetic variants from NGS reads.


[scanpy](scanpy.html) (1.8.1) Scanpy is a scalable toolkit for analyzing single-cell gene expression data. It includes preprocessing, visualization, clustering, pseudotime and trajectory inference and differential expression testing. The Python-based implementation efficiently deals with datasets of more than one million cells.


[schicexplorer](https://schicexplorer.readthedocs.io/) (7) scHiCExplorer is a set of programs to process, normalize, analyse and visualize single-cell Hi-C data


[Scramble](Scramble.html) (1.0.2) Scramble is a mobile element insertion (MEI) detection tool. It identifies clusters of soft clipped reads
in a BAM file, builds consensus sequences, aligns to representative L1Ta, AluYa5, and SVA-E sequences, and
outputs MEI calls.


[seqan](seqan.html) (2.4.0)  SeqAn is an open source C++ library of efficient algorithms and data structures for the analysis of sequences with the focus on biological data. It applies a unique generic design that guarantees high performance, generality, extensibility, and integration with other libraries. This package also contains a suite of apps, including Fiona, Gustaf, Mason, RazerS 3, Yara, SeqAn T-Coffee, Stellar, and searchjoin.


[seqkit](seqkit.html) (2.2.0) A cross-platform toolkit for FASTA/Q file manipulation


[seqtk](seqtk.html) (1.4) seqtk is a toolkit for processing sequences in FASTA/Q formats


[sequenza-utils](sequenza-utils.html) (3.0.0) Sequenza-utils is The supporting python library for the sequenza R package.

Sequenza is a project the estimate purity/ploidy and copy number alteration from tumor sequencing experiments. Sequenza-utils provide command lines programs to transform common NGS file type, such as BAM, pileup and VCF, to input files for the R package




[shapeit](shapeit.html) (5.1.0) SHAPEIT is a fast and accurate haplotype inference software


[shasta](shasta.html) (0.11.1) De novo assembly from Oxford Nanopore reads


[shmlast](shmlast.html) (1.6) shmlast is a reimplementation of the Conditional Reciprocal Best Hits algorithm for finding potential orthologs between a transcriptome and a species-specific protein database. It uses the LAST aligner and the pydata stack to achieve much better performance while staying in the Python ecosystem.



[shrimp](shrimp.html) (2\_2\_3) SHRiMP is a software package for aligning genomic reads against a target genome. It was primarily developed with the multitudinous short reads of next generation sequencing machines in mind, as well as Applied Biosystem's colourspace genomic representation. 


[sicer](sicer.html) (2-1.0.3) A clustering approach for identification of enriched domains from histone modification ChIP-Seq data


[sickle](sickle.html) (1.33) A windowed adaptive trimming tool for FASTQ files using quality


[SIFT](SIFT.html) (6.2.1) SIFT predicts whether an amino acid substitution affects protein function based on sequence homology and the physical properties of amino acids.


[slamdunk](slamdunk.html) (0.4.3) SlamDunk is a novel, fully automated software tool for automated, robust, scalable and reproducible SLAMseq data analysis. 


[smrtanalysis](smrtanalysis.html) (12.0.0.177059) SMRT® Analysis is a bioinformatics software suite available for analysis of DNA sequencing data from Pacific Biosciences’ SMRT technology. Users can choose from a variety of analysis protocols that utilize PacBio® and third-party tools. Analysis protocols include de novo genome assembly, cDNA mapping, DNA base-modification detection, and long-amplicon analysis to determine phased consensus sequences. 


[sniffles](sniffles.html) (2.2) Sniffles is a structural variation caller using third generation sequencing (PacBio or Oxford Nanopore). It detects all types of SVs (10bp+) using evidence from split-read alignments, high-mismatch regions, and coverage analysis.


[snippy](https://hpc.nih.gov/apps/snippy.html) (4.6.0) snippy is a tool for rapid haploid variant calling and core genome alignment


[snp2hla](snp2hla.html) (1.0.3) SNP2HLA is a tool to impute amino acid polymorphisms and single nucleotide polymorphisms in human luekocyte antigenes (HLA) within the major histocompatibility complex (MHC) region in chromosome 6. 


[snpEff](snpEff.html) (5.1d) snpEff is a variant annotation and effect prediction tool. It annotates and predicts the effects of variants on genes (such as amino acid changes).


snptest (2.5.6) SNPTEST is a program for the analysis of single SNP association in genome-wide studies. The tests implemented include

 \* Binary (case-control) phenotypes, single and multiple quantitative phenotypes
 \* Bayesian and Frequentist tests
 \* Ability to condition upon an arbitrary set of covariates
 \* Various different methods for the dealing with imputed SNPs.

The program is designed to work seamlessly with the output of both the genotype calling program CHIAMO, the genotype imputation program IMPUTE and the program GTOOL. 


[somaticsniper](somaticsniper.html) (1.0.5.0) The purpose of this program is to identify single nucleotide positions that are different between tumor and normal (or, in theory, any two bam files). It takes a tumor bam and a normal bam and compares the two to determine the differences. It outputs a file in a format very similar to Samtools consensus format. 


[sortmeRNA](sortmeRNA.html) (4.3.6) SortMeRNA is a local sequence alignment tool for filtering, mapping and clustering.
The core algorithm is based on approximate seeds and allows for sensitive analysis
of NGS reads. The main application of SortMeRNA is filtering rRNA
from metatranscriptomic data. SortMeRNA takes as input a file of reads
(fasta or fastq format) and one or multiple rRNA database file(s),
and sorts apart aligned and rejected reads into two files specified by the user.


[spaceranger](spaceranger.html) (2.1.1) 10x pipeline for processing Visium spatial RNA-seq data


[spades](spades.html) (3.15.5) SPAdes – St. Petersburg genome assembler – is intended for both standard isolates and single-cell MDA bacteria assemblies. 


[speedseq](speedseq.html) (0.1.2-20180208-4e60002) SpeedSeq is a genome analysis platform designed for rapid whole-genome variant detection and interpretation 


[sratoolkit](sratoolkit.html) (3.0.2) The NCBI SRA Toolkit enables reading ("dumping") of sequencing files from the SRA database and writing ("loading") files into the .sra format.


[SRST2](https://hpc.nih.gov/apps/srst2.html) (0.2.0) SRST2 is a a read mapping-based tool for rapid molecular typing of bacterial pathogens.
It allows fast and accurate detection of genes, alleles and multi-locus
sequence types (MLST) from WGS data. SRST2 is highly accurate and outperforms assembly-based methods in terms of both gene detection and allele assignment.


[stampy](stampy.html) (1.0.32) Short read aligner


[STAR](STAR.html) (2.7.10b) Spliced Transcripts Alignment to a Reference


[STAR-Fusion](STAR.html) (1.12.0) Transcript fusion detection


[STREAM](STREAM.html) (20180816) STREAM stands for Single-cell Trajectories Reconstruction, Exploration And Mapping ofomics data. It is an interactive pipeline capable of disentangling and visualizing
complex branching trajectories from both single-cell transcriptomic and epigenomic data.


[strelka](strelka.html) (2.9.10) Strelka is an analysis package designed to detect somatic SNVs and small indels from the aligned sequencing reads of matched tumor-normal samples.


[stringtie](stringtie.html) (2.2.1) StringTie is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts. It is primarily a genome-guided transcriptome assembler, although it can borrow algorithmic techniques from de novo genome assembly to help with transcript assembly. 


[subread](subread.html) (2.0.3) High-performance read alignment, quantification and mutation discovery


suppa (2.3) Fast, accurate, and uncertainty-aware differential splicing analysis across multiple conditions 


[SURVIVOR](survivor.html) (1.0.7) SURVIVOR is a tool set for simulating/evaluating SVs, merging and comparing SVs within and among samples, and includes various methods to reformat or summarize SVs.


[svanna](svanna.html) (1.0.3) The svanna is an efficient and accurate pathogenicity prediction for coding and regulatory structural variants in long-read genome sequencing. 


[svtyper](svtyper.html) (0.7.1) Svtyper is a Bayesian genotyper for structural variants.


[svviz](svviz.html) (1.6.2) svviz visualizes high-throughput sequencing data relevant to a structural variant. Only reads supporting the variant or the reference allele will be shown. svviz can operate in both an interactive web browser view to closely inspect individual variants, or in batch mode, allowing multiple variants (annotated in a VCF file) to be analyzed simultaneously. 


[talon](talon.html) (6.0) TALON is a Python package for identifying and quantifying known and novel genes/isoforms in long-read transcriptome data sets.


[taxonkit](taxonkit.html) (0.12.0) A Cross-platform and Efficient NCBI Taxonomy Toolkit



[telescope](telescope.html) (6cd5525) Single locus resolution of Transposable ELEment expression. Telescope estimates transposable element expression (retrotranscriptome) resolved to specific genomic locations. It directly addresses uncertainty in fragment assignment by reassigning ambiguously mapped fragments to the most probable source transcript as determined within a Bayesian statistical model.


telseq (0.0.2) TelSeq is a software that estimates telomere length from whole genome sequencing data (BAMs).


[tetoolkit](tetoolkit.html) (2.2.1) A package for including transposable elements in differential enrichment analysis of sequencing datasets. 


[THetA](theta.html) (0.7-20-g94fd772) Tumor Heterogeneity Analysis (THetA) is an algorithm used to estimate tumor purity and clonal/subclonal copy number aberrations simultaneously from high-throughput DNA sequencing data.


[tombo](tombo.html) (1.5.1) a suite of tools primarily for the identification of modified nucleotides from nanopore sequencing data. 


[tophat](tophat.html) (2.1.2) TopHat is a fast splice junction mapper for RNA-Seq reads. It aligns RNA-Seq reads to mammalian-sized genomes using the ultra high-throughput short read aligner Bowtie, and then analyzes the mapping results to identify splice junctions between exons. 


[TPMCalculator](TPMCalculator.html) (0.0.4) TPMCalculator quantifies mRNA abundance directly from the alignments by parsing BAM files


[transcript\_clean](transcript_clean.html) (2.0.3) TranscriptClean is a Python program that corrects mismatches, microindels, and noncanonical splice junctions in long reads that have been mapped to the genome. It is designed for use with sam files from the PacBio Iso-seq and Oxford Nanopore transcriptome sequencing technologies


[TransDecoder](TransDecoder.html) (5.5.0) TransDecoder identifies candidate coding regions within transcript sequences, such as those generated by de novo RNA-Seq transcript assembly using Trinity, or constructed based on RNA-Seq alignments to the genome using Tophat and Cufflinks.


[transvar](transvar.html) (2.5.9) TransVar is a versatile annotator for 3-way conversion and annotation among genomic characterization(s) of mutations and transcript-dependent annotation(s). 


[trimAl](trimAl.html) (1.2rev59) trimAl is a tool for the automated removal of spurious sequences or poorly aligned regions from a multiple sequence alignment.
It can consider several parameters, alone or in multiple combinations,
in order to select the most-reliable positions in the alignment.
These include the proportion of sequences with a gap, the level of residue similarity and,
if several alignments for the same set of sequences are provided, the consistency level of columns among alignments.
Moreover, trimAl is able to manually select a set of columns to be removed from the alignment.


[trimgalore](trimgalore.html) (0.6.7) Consistent quality and adapter trimming for RRBS or standard FastQ files.


[trimmomatic](trimmomatic.html) (0.39) Trimmomatic performs a variety of useful trimming tasks for illumina paired-end and single ended data.


[trinity](trinity.html) (2.15.1) Trinity, developed at the Broad Institute and the Hebrew University of Jerusalem, represents a novel method for the efficient and robust de novo reconstruction of transcriptomes from RNA-seq data.


[tRNAscan-SE](tRNAscan-SE.html) (2.0.9) tRNAscan-SE 2.0 has advanced the state-of-the-art methodology in tRNA gene detection and functional prediction, captured by rich new content of the companion Genomic tRNA Database


[ultraplex](ultraplex.html) (1.2.5) Ultraplex is primarily designed for the demultiplexing of sequencing data generated using in-house library preparation protocols with custom adaptors


[umitools](umitools.html) (1.1.2) tools for dealing with Unique Molecular Identifiers (UMIs)/Random Molecular Tags (RMTs) and single cell RNA-Seq cell barcodes


[vcf2maf](vcf2maf.html) (1.6.21) A smarter, more reproducible, and more configurable tool for converting a VCF to a MAF.


[vcfanno](vcfanno.html) (0.3.3) annotate a VCF with other VCFs/BEDs/tabixed files


[vcflib](vcflib.html) (1.0.3) a simple C++ library for parsing and manipulating VCF files, + many command-line utilities


[vcftools](vcftools.html) (0.1.16) VCFtools contains a Perl API (Vcf.pm) and a number of Perl scripts that can be used to perform common tasks with VCF files such as file validation, file merging, intersecting, complements, etc. 


[VEP](VEP.html) (110) VEP (Variant Effect Predictor) determines the effect of your variants (SNPs, insertions, deletions, CNVs or structural variants) on genes, transcripts, and protein sequence, as well as regulatory regions.


[verifybamid](verifybamid.html) (2.0.1) verifyBamID is a software that verifies whether the reads in particular file match previously known genotypes for an individual (or group of individuals), and checks whether the reads are contaminated as a mixture of two samples. verifyBamID can detect sample contamination and swaps when external genotypes are available. When external genotypes are not available, verifyBamID still robustly detects sample swaps. 


[verkko](verkko.html) (1.4.1) Verkko is a hybrid genome assembly pipeline developed for telomere-to-telomere assembly of PacBio HiFi and Oxford Nanopore reads.


[VirSorter2](VirSorter2.html) (2.2.3) VirSorter2 is a DNA and RNA virus identification tool. It leverages genome-informed
database advances across a collection of customized automatic classifiers
to improve the accuracy and range of virus sequence detection.


[VIRTUS](VIRTUS.html) (2.0.1) Bioinformatics pipeline for viral transcriptome detection.


[vsearch](vsearch.html) (2.22.1) VSEARCH supports de novo and reference based chimera detection, clustering, full-length and prefix dereplication, rereplication, reverse complementation, masking, all-vs-all pairwise global alignment, exact and global alignment searching, shuffling, subsampling and sorting. It also supports FASTQ file analysis, filtering, conversion and merging of paired-end reads.


[vt](vt.html) (0.57721) vt is a variant tool set that discovers short variants from Next Generation Sequencing data.



[winnowmap](winnowmap.html) (2.03) winnowmap is used for mapping ONT and PacBio reads to repetitive reference sequences.


[xHLA](xHLA.html) (2018-04-04) The HLA gene complex on human chromosome 6 is one of the
most polymorphic regions in the human genome and contributes
in large part to the diversity of the immune system. Accurate typing of HLA genes with short-read sequencing data has historically been difficult due to the sequence similarity between the polymorphic alleles. xHLA iteratively refines the mapping results at the amino acid level to
achieve high typing accuracy for both class I and II HLA genes.


[XHMM](XHMM.html) (2016-01-04) XHMM uses principal component analysis (PCA) normalization and a hidden Markov model (HMM) to detect and genotype copy number variation (CNV) from normalized read-depth data from targeted sequencing experiments.


[yak](yak.html) (r69) Yet another K-mer analyzer. Yak is initially developed for two specific use cases: 1) to robustly estimate the base accuracy of CCS reads and assembly contigs, and 2) to investigate the systematic error rate of CCS reads.


Image Analysis
[back to top](/apps/)  


[3DSlicer](3DSlicer.html) (5.2.2) A software platform for the analysis (including registration and interactive segmentation) and visualization (including volume rendering) of medical images and for research in image guided therapy. 


[AFNI](afni.html) (current-py3) AFNI (Analysis of Functional NeuroImages) is a set of C programs for processing, analyzing, and displaying functional MRI (FMRI) data - a technique for mapping human brain activity.




[ANTs](ANTs.html) (2.4.2) Advanced Normalization Tools (ANTs) extracts information from complex datasets that include imaging. Paired with ANTsR (answer), ANTs is useful for managing, interpreting and visualizing multidimensional data. 


AreTomo (1.3.4) Alignment and Reconstruction for Electron Tomography


[B-SOID](b-soid.html) (1.3) B-SOID (Behavioral Segmentation in Deeplabcut) is an 
unsupervised learning algorithm that serves to discover and classify behaviors that are not pre-defined by users.
It segregates statistically different, sub-second rodent behaviors with a single bottom-up perspective video cameraR by performing a novel expectation maximization
fitting of Gaussian mixture models on t-Distributed Stochastic Neighbor Embedding (t-SNE).


[baracus](https://hpc.nih.gov/apps/baracus.html) (1.1.4) Baracus predicts brain age, based on data from Freesurfer. It combines data from cortical thickness, cortical surface area, and subcortical information 


[BrkRaw](BrkRaw.html) (0.3.6) The ‘BrkRaw’ is a python module designed to provide a comprehensive tool to access raw data acquired from Bruker Biospin preclinical MRI scanner. This module is also compatible with the zip compressed data to enable use of the archived data directly.
The module is comprised of four components, including graphical user interface (GUI), command-line tools, high-level and low-level python APIs.


[Bsoft](Bsoft.html) (1.9.0) Bsoft is a collection of programs and a platform for development of software for image and molecular processing in structural biology. Problems in structural biology are approached with a highly modular design, allowing fast development of new algorithms without the burden of issues such as file I/O. It provides an easily accessible interface, a resource that can be and has been used in other packages.



[c3d](https://hpc.nih.gov/apps/c3d.html) (1.1.0) C3D is a command-line tool for converting 3D images between common file formats. The tool also includes a growing list of commands for image manipulation, such as thresholding and resampling.


[cellpose](cellpose.html) (2.1.1) A generalist algorithm for cellular segmentation with human-in-the-loop capabilities.


[cisTEM](cisTEM.html) (1.0.0-beta) cisTEM is user-friendly software to process cryo-EM images of macromolecular complexes and obtain high-resolution 3D reconstructions from them.


[civet](https://hpc.nih.gov/apps/civet.html) (2.1.1) civet is a brain-imaging pipeline for
analysis of large MR data sets. civet
extracts and analyses cortical surfaces 
from MR images, as well as many other 
volumetric and corticometric functions.


[cmtk](cmtk.html) (3.3.2) CMTK is a Software toolkit for computational morphometry of biomedical images. CMTK provides a set of command line tools for processing and I/O.


CompuCell3D (4.4.0) CompuCell3D is a multiscale multicellular virtual tissue modeling and simulation environment. CompuCell3D is written in C++ and provides Python bindings for model and simulation development in Python. CompuCell3D is supported on Windows, Mac and Linux.


[connectome-workbench](connectome-workbench.html) (1.5.0) Tools to browse, download, explore, and analyze data from the Human Connectome Project (HCP). Allows users to compare their own data to that of the HCP. 


[coolbox](coolbox.html) (0.3.8) CoolBox is an open-source, user-friendly toolkit for visual analysis of genomics data.
It is highly compatible with the Python ecosystem and customizable with a well-designed user interface.
It can bed used, for example, to produce high-quality genome track
plots or fetch commonly used genomic data files with a Python script or command
line, to explore genomic data interactively within Jupyter environment or web browser.


[cryosparc](https://hpc.nih.gov/nih/cryosparc/) (4.2.1) CryoSPARC (Cryo-EM Single Particle Ab-Initio Reconstruction and Classification) is a state of the art HPC software solution for complete processing of single-particle cryo-electron microscopy (cryo-EM) data. CryoSPARC is useful for solving cryo-EM structures of membrane proteins, viruses, complexes, flexible molecules, small particles, phase plate data and negative stain data.


[CTF](afni.html) (6.1) The CTF MEG software has two main roles:
- Provide a human-machine interface to the CTF MEG elec- tronics to collect MEG and/or EEG data.
- Provide a tool for reviewing and (to a limited extent) ana- lyzing the MEG and/or EEG data acquired by the CTF MEG system.


[ctffind](ctffind.html) (4.1.14) Programs for finding CTFs of electron micrographs


[dcm2niix](dcm2niix.html) (1.0.20211006) DICOM to NIfTI converter


[DeepCAD](DeepCAD.html) (20210826) DeepCAD is a self-supervised deep-learning method
for spatiotemporal enhancement of calcium imaging data
that does not require any high signal-to-noise ratio (SNR) observations.
DeepCAD suppresses detection noise and improves the SNR more than tenfold,
which reinforces the accuracy of neuron extraction and spike inference
and facilitates the functional analysis of neural circuits.


[DeepLabCut](https://hpc.nih.gov/apps/DeepLabCut.html) (2.3.4) DeepLabCut is an open source toolbox
that builds on a state-of-the-art human pose estimation algorithm.
It allows training of a deep neural network by using limited training data
to precisely track user-defined features, so that the human labeling accuracy will be matched.


[deepmedic](deepmedic.html) (0.8.4) This project aims to offer easy access to Deep Learning for segmentation of structures of interest in biomedical 3D scans. It is a system that allows the easy creation of a 3D Convolutional Neural Network, which can be trained to detect and segment structures if corresponding ground truth labels are provided for training. The system processes NIFTI images, making its use straightforward for many biomedical tasks.


[dmriprep](dmriprep.html) (0.5.0) Pipeline for preprocessing diffusion MRI datasets.


[dynamo](https://hpc.nih.gov/apps/dynamo.html) (1.1.532) Dynamo is a software environment for subtomogram averaging of cryo-EM data.


[elastix](elastix.html) (4.9; 5.1.0) a toolbox for rigid and nonrigid registration of images.


[EMAN2](EMAN2.html) (2.99) EMAN2 is a broadly based greyscale scientific image processing suite with a primary focus on processing data from transmission electron microscopes. 


[fastsurfer](https://hpc.nih.gov/apps/fastsurfer.html) (1.1.1) Fastsurfer is a neuroimaging pipeline based on deep learning.


[Fiji](Fiji.html) (1.52f) Fiji is a distribution of the popular open-source software ImageJ focused on biological-image analysis. Fiji uses modern software engineering practices to combine powerful software libraries with a broad range of scripting languages to enable rapid prototyping of image-processing algorithms. Fiji facilitates the transformation of new algorithms into ImageJ plugins that can be shared with end users through an integrated update system.


[fitlins](https://hpc.nih.gov/apps/fitlins.html) (0.9.1) Fitlins fits linear models to BIDS neuroimaging datasets.


[fmriprep](https://hpc.nih.gov/apps/fmriprep.html) (23.1.4) A Robust Preprocessing Pipeline for fMRI Data


[Frealign](Frealign.html) (9.11\_151031) Frealign is a program for high-resolution refinement of 3D reconstructions from cryo-EM images of single particles.


[Freesurfer](freesurfer.html) (7.4.1) Freesurfer is a set of automated tools for reconstruction of the brain's cortical surface from structural MRI data, and overlay of functional MRI data onto the reconstructed surface.


[FSL](fsl.html) (6.0.6) FSL is a comprehensive library of image analysis and statistical tools for FMRI, MRI and DTI brain imaging data.


[Gctf](Gctf.html) (1.06) Gctf provides accurate estimation of the contrast transfer function (CTF) for near-atomic resolution cryo electron microscopy (cryoEM) reconstruction using GPUs. 


[Huygens](Huygens.html) (23.04.0-p6) Huygens is an image restoration, deconvolution, resolution and noise reduction. It can process images from all current optical microscopes, including wide-field, confocal, Nipkow (scanning disk confocal), multiple-photon, and 4Pi microscopes.


[Imaris](imaris.html) (9.7) Imaris provides scientists with solutions for processing, visualizing and analyzing multi-dimensional microscopic images. It reads images in many of the most commonly used proprietary formats.


[IMOD](IMOD.html) (4.12.25) IMOD is a set of image processing, modeling and display programs used for tomographic reconstruction and for 3D reconstruction of EM serial sections and optical sections. 


[IsoNet](IsoNet.html) (0.2.1) IsoNet is a deep learning-based software package that iteratively reconstructs the missing-wedge information and increases signal-to-noise ratio, using the knowledge learned from raw tomograms. Without the need for sub-tomogram averaging, IsoNet generates tomograms with significantly reduced resolution anisotropy.


ITK-SNAP (3.8.0) ITK-SNAP is a tool for segmentation of 3D biomedical images. It requires a graphical connection to run on the cluster. 


[ITK-SNAP-BS](ITK-SNAP-BS.html) (20131007) ITK-SNAP brings active contour
segmentation to the fingertips of clinical researchers.
This application fulfills a specific and pressing need of biomedical
imaging research by providing a combination of manual and
semiautomatic tools for extracting structures in 3D image data of
different modalities and from different anatomical regions.


laynii (2.4.0) Tools to analyze layer fMRI datasets


magetbrain (1.0) Given a set of labelled MR images (atlases) and unlabelled images (subjects), MAGeT produces a segmentation for each subject using a multi-atlas voting procedure based on a template library made up of images from the subject set.


[mango](https://hpc.nih.gov/apps/mango.html) (4.1) Mango (Multi-image Analysis GUI) is a viewer for medical research images. It provides analysis tools and a user interface to navigate image volumes. 


[minc-toolkit](https://github.com/BIC-MNI/minc-toolkit-v2) (1.9.18) This metaproject bundles multiple MINC-based packages that historically have been developed somewhat independently


[MIPAV](MIPAV.html) (11.0.3) The MIPAV (Medical Image Processing, Analysis, and Visualization) application enables quantitative analysis and visualization of medical images of numerous modalities such as PET, MRI, CT, or microscopy. 


[MONAILabel](MONAILabel.html) (0.3.2; 0.4.2) MONAI Label is a free and open-source platform that facilitates the development of AI-based applications
that aim at reducing the time required to annotate 3D medical image datasets.
It allows researchers to readily deploy their apps as services, which can be made available to clinicians
via their preferred user-interface. Currently, MONAI Label readily supports locally installed (3DSlicer) and web-based (OHIF) frontends, and offers two Active learning strategies to facilitate and speed up the training of segmentation algorithms. MONAI Label allows researchers to make incremental improvements to their labeling apps by making them available to other researchers and clinicians alike.


[MotionCor2](MotionCor2.html) (1.5.0) MotionCor2 is a multi-GPU accelerated program that provides iterative, patch-based motion detection combining spatial and temporal constraints and dose weighting for both single particle and tomographic cryo-electon microscopy images.


mricron (1.0.20190902) Viewer for several types of brain scan formats


[mriqc](https://hpc.nih.gov/apps/mriqc.html) (23.1.0) MRIQC is an MRI quality control tool


[mrtrix](mrtrix.html) (3.0.4) MRtrix provides a large suite of tools for image processing, analysis and visualisation, with a focus on the analysis of white matter using diffusion-weighted MRI.


[mrtrix3tissue](mrtrix3tissue.html) (5.2.9) MRtrix3Tissue is a fork of the MRtrix3 project. It aims to add capabilities for 3-Tissue CSD modelling and analysis to a complete version of the MRtrix3 software.


napari (0.4.18) napari is an image viewer.Can be used to annotate images.


[nibabies](https://hpc.nih.gov/apps/nibabies.html) (22.0.2) Preprocessing pipeline for neonate and infant MRI.


OpenSlide (3.4.1) OpenSlide is a C library for reading and manipulating digital slides of diverse
vendor formats. It provides a simple interface to read whole-slide images
(also known as virtual slides). OpenSlide has been used in the digital pathology projects.


[PEET](PEET.html) (1.15.0) PEET (Particle Estimation for Electron Tomography) is an open-source package for aligning and averaging particles in 3-D subvolumes extracted from tomograms. It seeks the optimal alignment of each particle against a reference volume through several iterations. If PEET and IMOD are both installed, most PEET operations are available from the eTomo graphical user interface in IMOD. 


[petprep\_hmc](petprep_hmc.html) (0.06) Positron Emission Tomography (PET) is a state-of-the-art neuroimaging tool for quantification of the in vivo spatial distribution of specific molecules in the brain. It is affected by various kinds of patient movement during a scan.
The PETprep\_HMC application allows for correction of the PET results in the presence of head motions.


[plastimatch](https://hpc.nih.gov/apps/plastimatch.html) (1.9.3) Application for registration of medical images such as X-rays, CT, MRI and PET


[PyTom](PyTom.html) (1.0) PyTom is a software package for the analysis of volumetric data
obtained by cryo electron tomography (cryo-ET).
It covers a complete pipeline of processing steps
for tomogram reconstruction, localization of macromolecular complexes in tomograms,
fine alignment of subtomograms extracted at these locations, and their classification.


[qsiprep](qsiprep.html) (0.16.1) qsiprep configures pipelines for processing diffusion-weighted MRI (dMRI) data.


[QuPath](QuPath.html) (0.4.1) QuPath is open source software for bioimage analysis. It is often used for digital pathology applications because it offers a powerful set of tools for working with whole slide images - but it can be applied to lots of other kinds of image as well.


[RELION](RELION/index.html) (4.0.1) RELION (for REgularised LIkelihood OptimisatioN) is a stand-alone computer program for Maximum A Posteriori refinement of (multiple) 3D reconstructions or 2D class averages in cryo-electron microscopy.


[ResMap](ResMap.html) (1.9) ResMap (Resolution Map) is a Python (NumPy/SciPy) application with a Tkinter GUI and a command-line interface. It is a software package for computing the local resolution of 3D density maps studied in structural biology, primarily electron cryo-microscopy (cryo-EM).


[SAMsrcV3](SAMsrcV3.html) (20180713-c5e1042) Synthetic Aperture Magnetometry - The SANsrcV3 suite implements the latest advances in MEG source localization. 


[SimNIBS](simnibs.html) (4.0) SimNIBS is a free software package for the Simulation of Non-invasive Brain Stimulation. It allows for realistic calculations of the electric field induced by transcranial magnetic stimulation (TMS) and transcranial direct current stimulation (tDCS).


[smriprep](smriprep.html) (0.8.3) Structural MRI PREProcessing (sMRIPrep) workflows for NIPreps (NeuroImaging PREProcessing tools).


[SPHIRE](SPHIRE.html) (1.4) SPHIRE (SPARX for High-Resolution Electron Microscopy) is an open-source,
user-friendly software suite for the semi-automated
processing of single particle electron cryo-microscopy (cryo-EM) data.
It allows fast and reproducible structure determination from cryo-EM images.


[spm12](https://hpc.nih.gov/apps/spm12.html) (7870) The (S)tatistical (P)ara(M)etric application analyzes brain imaging data.


tedana (23.0.1) Tedana is an application to denoise multi-echo fMRI datasets


[tomotwin](tomotwin.html) (0.5.0) TomoTwin - a deep metric learning based particle picking procedure for cryo-ET


[topaz](https://hpc.nih.gov/apps/topaz.html) (0.2.5-ba91e19) topaz is a pipeline for particle detection in 
cryo-electron microscopy images using convolutional 
neural networks trained from positive and unlabeled 
examples. Topaz also includes methods for micrograph 
and tomogram denoising using deep denoising models.



[TORTOISE](TORTOISE.html) (3.2.0) (Tolerably Obsessive Registration and Tensor Optimization
Indolent Software Ensemble)
The TORTOISE software package is for processing diffusion MRI data.


[tractseg](tractseg.html) (1.7.1) Tool for white matter bundle segmentation from Diffusion MRI.


[xcpengine](https://hpc.nih.gov/apps/xcpengine.html) (1.2.4) xpcEngine performs denoising and estimation of Functional Connectivity on fMRI datasets


[xcp\_d](xcp_d.html) (0.5.0) xcp\_d is a postprocessing and noise regression pipeline for fMRI datasets (can use output from fmriprep and nibabies).


Linkage/Phylogenetics
[back to top](/apps/)  


[AdmixTools](AdmixTools.html) (7.0.2) ADMIXTOOLS is a software package that supports formal tests of whether admixture occurred, and makes it possible to infer admixture proportions and dates. 


[admixture](admixture.html) (1.3.0) ADMIXTURE is a software tool for maximum likelihood estimation of individual ancestries from multilocus SNP genotype datasets. It uses the same statistical model as STRUCTURE but calculates estimates much more rapidly using a fast numerical optimization algorithm. 


[arcashla](arcashla.html) (0.5.0) high resolution HLA typing from RNA seq


[Beagle](Beagle.html) (5.4\_22Jul22) Beagle is a package for imputing genotypes, inferring haplotype phase, and performing genetic association analysis. BEAGLE is designed to analyze large-scale data sets with hundreds of thousands of markers genotyped on thousands of samples.


[beagle-lib](https://beagle-dev.github.io/html/beagle_8h.html) (4.0) BEAGLE is a high-performance library that can perform the core calculations at the heart of most Bayesian and Maximum Likelihood phylogenetics packages. It can make use of highly-parallel processors such as those in graphics cards (GPUs) found in many PCs.
  

module name: beagle-lib


[BEAST](BEAST.html) (1.10.4,2.6.2) BEAST (Bayesian Evolutionary Analysis Sampling Trees) is a cross-platform program for Bayesian MCMC analysis of molecular sequences.


[biobakery\_workflows](biobakery_workflows.html) (3.1) bioBakery is a meta’omic analysis environment and collection of individual software
tools with the capacity to process raw shotgun sequencing data into actionable microbial community
feature profiles, summary reports, and publication-ready figures.
It includes a collection of preconfigured analysis modules also joined into workflows for reproducibility.
Each individual module has been developed to perform a particular task, e.g. quantitative taxonomic
profiling or statistical analysis.


[CD-HIT](cd-hit.html) (4.6.8) CD-HIT is a very widely used program for clustering and comparing protein or nucleotide sequences.


[eigensoft](eigensoft.html) (8.0.0) The EIGENSOFT package combines functionality from population genetics methods and EIGENSTRAT stratification correction method.


[FastQTL](FastQTL.html) (2.184) In order to discover quantitative trait loci (QTLs),
multi-dimensional genomic datasets combining
DNA-seq and ChiP-/RNA-seq require methods that rapidly correlate tens of thousands of
molecular phenotypes with millions of genetic variants while appropriately controlling for multiple testing. FastQTL implements a popular cis-QTL mapping strategy
in a user- and cluster-friendly tool. FastQTL also proposes an efficient
permutation procedure to control for multiple testing.



[FastTree](FastTree.html) (2.1.11) FastTree infers approximately-maximum-likelihood phylogenetic trees from alignments of nucleotide or protein sequences. FastTree can handle alignments with up to a million of sequences in a reasonable amount of time and memory. For large alignments, FastTree is 100-1,000 times faster than PhyML 3.0 or RAxML 7. 


[GCTA](GCTA.html) (1.94.0beta0) GCTA (Genome-wide Complex Trait Analysis) is designed to estimate the proportion of phenotypic variance explained by genome- or chromosome-wide SNPs for complex traits.


[gtdb-tk](gtdb-tk.html) (2.3.2) GTDB-Tk is a software toolkit for assigning objective taxonomic classifications to bacterial and archaeal genomes based on the Genome Database Taxonomy


[gubbins](https://hpc.nih.gov/apps/gubbins.html) (2.3.4) Gubbins is an algorithm that iteratively identifies loci containing elevated densities of base substitions while concurrently constructing a phylogeny based on the putative point mutations outside of these regions.



[hyphy](hyphy.html) (2.5.29) HyPhy (Hypothesis Testing using Phylogenies) is an open-source software package for the analysis of genetic sequences (in particular the inference of natural selection) using techniques in phylogenetics, molecular evolution, and machine learning.


[iqtree](iqtree.html) (2.2.0.5) Efficient phylogenomic software by maximum likelihood


[king](king.html) (2.2.7) KING is a toolset to explore genotype data from a genome-wide association study (GWAS) or a sequencing project. KING can be used to check family relationship and flag pedigree errors by estimating kinship coefficients and inferring IBD segments for all pairwise relationships.


[LTSOFT](LTSOFT.html) (4.0) The LTSOFT application implements a new approach
to using information from known associated variants
when conducting disease association studies.
The approach is based in part on the classical technique
of liability threshold modeling and performs estimation
of model parameters for each known variant
while accounting for the published disease prevalence from
the epidemiological literature.


[Madeline2](Madeline.html) (2.0) The Madeline 2.0 Pedigree Drawing Engine (PDE) is a
pedigree drawing program for use in linkage and family-based
association studies. The program is designed to handle large and
complex pedigrees with an emphasis on readability and aesthetics.


[MEGA](MEGA.html) (11.0.10) MEGA, Molecular Evolutionary Genetics Analysis, is a software suite for analyzing DNA and protein sequence data from species and populations


[mega2](mega2.html) (6.0.0) Mega2 is a data-handling program for facilitating genetic linkage and association analyses.


[merlin](merlin.html) (1.1.2) MERLIN uses sparse trees to represent gene flow in pedigrees and is one of the fastest pedigree analysis packages around.


[mothur](mothur.html) (1.48.0) mothur is a tool for analyzing 16S rRNA gene sequences generated on multiple platforms as part of microbial ecology projects.


[MrBayes](mrbayes.html) (3.2.7) MrBayes is a program for Bayesian inference and model choice across a wide range of phylogenetic and evolutionary models. MrBayes uses Markov chain Monte Carlo (MCMC) methods to estimate the posterior distribution of model parameters.


[ohana](ohana.html) (current) Ohana is a suite of software for analyzing population structure and admixture history using unsupervised learning methods.


[phase](phase.html) (2.1.1) infers haplotypes from population genotype data


[Phylip](Phylip.html) (3.698) Phylip is a package of programs for inferring phylogenies (evolutionary trees). Includes methods for parsimony, distance matrix and likelihood methods.




[pplacer](pplacer.html) (1.1) Pplacer places query sequences on a fixed reference phylogenetic tree to maximize phylogenetic likelihood or posterior probability according to a reference alignment. Pplacer is designed to be fast, to give useful information about uncertainty, and to offer advanced visualization and downstream analysis. 


[PRIMUS](https://hpc.nih.gov/apps/PRIMUS.html) (1.9.0) PRIMUS, Pedigree Reconstruction and Identification of a Maximum Unrelated Set, is used to read genome-wide IBD estimates and identify a maximum unrelated set and pedigree reconstruction.


[PRSice](prsice.html) (2.3.3) PRSice is a Polygenic Risk Score software for calculating, applying, evaluating and plotting the results of polygenic risk scores (PRS) analyses. 


[QIIME](QIIME.html) (2023.5) QIIME is an open source software package for comparison and analysis of microbial communities, primarily based on high-throughput amplicon sequencing data (such as SSU rRNA) generated on a variety of platforms, but also supporting analysis of other types of data (such as shotgun metagenomic data). 


[RAxML](raxml.html) (8.2.12) RAxML (randomized axelerated maximum likelihood) is a sequential and parallel program for inference of large phylogenies with maximum likelihood (ML). 


[raxml-ng](raxml-ng.html) (1.2.0) RAxML-NG is a phylogenetic tree inference tool which uses maximum-likelihood (ML) optimality criterion. Its search heuristic is based on iteratively performing a series of Subtree Pruning and Regrafting (SPR) moves, which allows to quickly navigate to the best-known ML tree. Successor to raxml.


[RevBayes](revbayes.html) (1.2.1) Bayesian phylogenetic inference using probabilistic graphical models and an interpreted language


[scDRS](scDRS.html) (1.02) The scDRS application implements an approach that links scRNA-seq
with polygenic disease risk at single-cell resolution, independent of annotated cell types.
scDRS identifies cells exhibiting excess expression across disease-associated genes
implicated by genome-wide association studies (GWASs).
Genes whose expression was correlated with the scDRS score across cells
(reflecting coexpression with GWAS disease-associated genes)
were strongly enriched for gold-standard drug target and Mendelian disease genes.


[shapeit](shapeit.html) (5.1.0) SHAPEIT is a fast and accurate haplotype inference software


[SMR](SMR.html) (0.702; 1.03; 1.3.1) SMR integrates summary-level data from GWAS with data from expression quantitative trait locus (eQTL) studies to identify genes whose expression levels are associated with a complex trait because of pleiotropy. It implements methods to test for pleiotropic association between the expression level of a gene and a complex trait of interest using summary-level data from GWAS and expression quantitative trait loci (eQTL) studies (Zhu et al. 2016 Nat Genet).


[solar](solar.html) (9.0.1) SOLAR-Eclipse is an extensive, flexible software package for genetic variance components analysis, including linkage analysis, quantitative genetic analysis, SNP association analysis (QTN and QTLD), and covariate screening.




[sumtrees](sumtrees.html) (4.5.2) Summarize non-parameteric bootstrap or Bayesian posterior probability support for splits or clades on phylogenetic trees.


[TensorQTL](TensorQTL.html) (1.0.7) ensoorQTL leverages general-purpose libraries and graphics processing units
(GPUs) to achieve high efficiency of computations at low costR. Using PyTorch or TensorFlow
it allows > 200-fold decreases in runtime and ~ 5–10-fold reductions in cost when
running on GPUs relative to CPUs.


[treemix](treemix.html) (1.13) TreeMix is a method for inferring the patterns of population splits and mixtures in the history of a set of populations.



[treetime](treetime.html) (0.11.1) TreeTime provides routines for ancestral sequence reconstruction and inference of molecular-clock phylogenies.


Mass Spectrometry
[back to top](/apps/)  
