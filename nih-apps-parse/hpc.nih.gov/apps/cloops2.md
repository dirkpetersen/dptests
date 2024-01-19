

document.querySelector('title').textContent = "cloops2";
CLOOPS2 on Biowulf


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



cLoops2: full stack analysis tool for chromatin interactions



### References:


* Yaqiang Cao, Shuai Liu, Gang Ren, Qingsong Tang, Keji Zhao
 [**cLoops2: a full-stack comprehensive analytical tool for chromatin interactions**](https://doi.org/10.1093/nar/gkab1233)
*Nucleic Acids Research, Volume 50, Issue 1, 11 January 2022, Pages 57–71*


Documentation
* [cLoops2 Main Site](https://github.com/YaqiangCao/cLoops2)


Important Notes
* Module Name: chromeister (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded
* Environment variables set 
	+ CLOOPS2\_HOME



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
[user@cn4274 ~]$ **module load cloops2**
[+] Loading cloops2  0.0.4  on cn4326
[+] Loading singularity  3.10.5  on cn4326

```


Example
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Execute cLoops2



```

[user@cn4274 ~]$ **cLoops2 -h**
An enhanced, accurate and flexible peak/domain/loop-calling and analysis tool
for 3D genomic interaction data.

Use cLoops2 sub-command -h to see detail options and examples for sub-commands.
Available sub-commands are:
    qc: quality control of BEDPE files before analysis.
    pre: preprocess input BEDPE files into cLoops2 data.
    update: update cLoops2 data files locations.
    combine: combine multiple cLooops2 data directories.
    dump: convert cLoops2 data files to others (BEDPE, HIC, washU, bedGraph,
          contact matrix or virtual 4C signal)
    estEps: estimate eps using Gaussian mixture models or k-distance plot.
    estRes: estimate reasonable contact matrix resolution based on signal
            enrichment.
    estDis: estimate significant interactions distance range.
    estSat: estimate sequencing saturation based on contact matrix.
    estSim: estimate similarities among samples based on contact matrix.
    filterPETs: filter PETs based on peaks, loops, singleton mode or knn mode.
    samplePETs: sample PETs according to specific target size.
    callPeaks: call peaks for ChIP-seq, ATAC-seq, ChIC-seq and CUT&Tag or the
               3D genomic data such as Trac-looping, Hi-TrAC, HiChIP and more.
    callLoops: call loops for 3D genomic data.
    callDiffLoops: call differentially enriched loops for two datasets.
    callDomains: call domains for 3D genomic data.
    plot: plot the interaction matrix, genes, view point plot, 1D tracks,
          peaks, loops and domains for a specific region.
    montage: analysis of specific regions, producing Westworld Season 3 -like
             Rehoboam plot.
    agg: aggregated feature analysis and plots, features can be peaks, view
         points, loops and domains.
    quant: quantify peaks, loops and domains.
    anaLoops: anotate loops for target genes.
    findTargets: find target genes of genomic regions through networks from
                 anaLoops.

Examples:
    cLoops2 qc -f trac_rep1.bedpe.gz,trac_rep2.bedpe,trac_rep3.bedpe.gz \
               -o trac_stat -p 3
    cLoops2 pre -f ../test_GM12878_chr21_trac.bedpe -o trac
    cLoops2 update -d ./trac
    cLoops2 combine -ds ./trac1,./trac2,./trac3 -o trac_combined -keep 1
    cLoops2 dump -d ./trac -o trac -hic
    cLoops2 estEps -d trac -o trac_estEps_gmm -p 10 -method gmm
    cLoops2 estRes -d trac -o trac_estRes -p 10 -bs 25000,5000,1000,200
    cLoops2 estDis -d trac -o trac -plot -bs 1000
    cLoops2 estSim -ds Trac1,Trac2 -o trac_sim -p 10 -bs 2000 -m pcc -plot
    cLoops2 filterPETs -d trac -peaks trac_peaks.bed -o trac_peaksFiltered -p 10
    cLoops2 samplePETs -d trac -o trac_sampled -t 5000000 -p 10
    cLoops2 callPeaks -d H3K4me3_ChIC -bgd IgG_ChIC -o H3K4me3_cLoops2 -eps 150 \
                      -minPts 10
    cLoops2 callLoops -d Trac -eps 200,500,1000 -minPts 3 -filter -o Trac -w -j \
                      -cut 2000
    cLoops2 callLoops -d HiC -eps 1000,5000,10000 -minPts 10,20,50,100 -w -j \
                      -trans -o HiC_trans
    cLoops2 callDiffLoops -tloop target_loop.txt -cloop control_loop.txt \
                          -td ./target -cd ./control -o target_diff
    cLoops2 callDomains -d trac -o trac -bs 10000 -ws 200000
    cLoops2 plot -f test/chr21-chr21.ixy -o test -bs 500 -start 34840000 \
                 -end 34895000 -triu -1D -loop test_loops.txt -log \
                 -gtf hg38.gtf -bws ctcf.bw -beds enhancer.bed
    cLoops2 montage -f test/chr21-chr21.ixy -o test -bed test.bed
    cLoops2 agg -d trac -loops trac.loop -peaks trac_peaks.bed \
                -domains hic_domains.bed -bws CTCF.bw,ATAC.bw -p 20 -o trac
    cLoops2 quant -d trac -peaks trac_peaks.bed -loops trac.loop \
                  -domains trac_domain.txt -p 20 -o trac
    cLoops2 anaLoops -loops test_loop.txt -gtf gene.gtf -net -o test
    cLoops2 findTargets -net test_ep_net.sif -tg test_targets.txt \
                        -bed GWAS.bed -o test
    More usages and examples are shown when run with cLoops2 sub-command -h.


optional arguments:
  -h, --help  show this help message and exit
  -d PREDIR   Assign data directory generated by cLoops2 pre to carry out analysis.
  -o FNOUT    Output data directory / file name prefix, default is cLoops2_output.
  -p CPU      CPUs used to run the job, default is 1, set -1 to use all CPUs
              available. Too many CPU could cause out-of-memory problem if there are
              too many PETs.
  -cut CUT    Distance cutoff to filter cis PETs, only keep PETs with distance
              >=cut. Default is 0, no filtering.
  -mcut MCUT  Keep the PETs with distance <=mcut. Default is -1, no filtering.
  -v          Show cLoops2 verison number and exit.
  ---         Following are sub-commands specific options. This option just show
              version of cLoops2.

Bug reports are welcome and can be put as issue at github repo or sent to
caoyaqiang0410@gmail.com or yaqiang.cao@nih.gov. Thank you.

```

Annotate an SV:



```

[user@cn4338] **cp -a /usr/local/apps/duphold/0.2.3/test\_data .**
[user@cn4338 test_data]$ **duphold \
 --threads 4 \
 --vcf sparse\_in.vcf \
 --bam sparse.cram \
 --fasta sparse.fa \
 --output output.bcf** 
#To view output, load samtools and view with bcftools
[user@cn4338 test_data] **module load samtools**
[user@cn4338 test_data] **bcftools view test-out.bcf**
##fileformat=VCFv4.2
...
##bcftools_viewVersion=1.4-19-g1802ff3+htslib-1.4-29-g42bfe70
##bcftools_viewCommand=view CHM1_CHM13/full.37d5.vcf.gz; Date=Mon Sep 24 13:48:04 2018
...
##bcftools_viewVersion=1.17+htslib-1.17
##bcftools_viewCommand=view test-out.bcf; Date=Thu May 25 12:49:34 2023
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Eluc-CR2.F
NW_017858824.1  135118  72454   N       DEL   5875.46 .       SVTYPE=DEL;END=135332;CIPOS=0,0;CIEND=0,0;CIPOS95=0,0;CIEND95=0,0;GCF=0.306977  GT:DP:DHFC:DHFFC:DHBFC:DHSP     0/1:200:1.91667:0.597403:1.76923:0

```



|  |  |
| --- | --- |
|  For more information on pre and post processing, please visit the [cLoops2 Github](https://github.com/YaqiangCao/cLoops2#-cloops2-full-stack-analysis-tool-for-chromatin-interactions-) | |








