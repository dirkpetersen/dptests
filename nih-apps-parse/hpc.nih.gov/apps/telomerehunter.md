

document.querySelector('title').textContent = 'TelomereHunter: estimation of telomere content and composition from cancer genomes. ';
**TelomereHunter: estimation of telomere content and composition from cancer genomes.** 


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



TelomereHunter is a software for the detailed characterization of telomere
maintenance mechanism footprints in the genome. The tool is implemented for the analysis of large cancer
genome cohorts and provides a variety of diagnostic diagrams as well as machine-readable output for
subsequent analysis. 



### References:


* Lars Feuerbach, Lina Sieverling, Katharina I. Deeg, Philip Ginsbach, Barbara Hutter, Ivo Buchhalter,
Paul A. Northcott, Sadaf S. Mughal, Priya Chudasama, Hanno Glimm, Claudia Scholl, Peter Lichter,
Stefan Fröhling, Stefan M. Pfister, David T. W. Jones, Karsten Rippe and Benedikt Brors  

 *TelomereHunter – in silico estimation of telomere content and composition from cancer genomes*   

 [BMC Bioinformatics](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2851-0) (2019) 20:272 https://doi.org/10.1186/s12859-019-2851-0.


Documentation
* [TelomereHunter Home page](https://www.dkfz.de/en/applied-bioinformatics/telomerehunter/telomerehunter.html)


Important Notes
* Module Name: eggNOG-mapper (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **TH\_HOME**  installation directory
	+ **TH\_BIN**       executable directory
	+ **TH\_SRC**       source code directory
	+ **TH\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=16g --gres=lscratch:30 -c 8**
[user@cn3200 ~]$**module load telomerehunter**
[user@cn3200 ~]$**telomerehunter -h**


        TelomereHunter  Copyright (C) 2015  Lina Sieverling, Philip Ginsbach, Lars Feuerbach
        This program comes with ABSOLUTELY NO WARRANTY.
        This is free software, and you are welcome to redistribute it
        under certain conditions. For details see the GNU General Public License
        in the license copy received with TelomereHunter or .


TelomereHunter 1.1.0


usage: telomerehunter [-h] [-ibt TUMOR_BAM] [-ibc CONTROL_BAM] -o OUTPUT_DIR
                      -p PID [-b BANDING_FILE] [-rt REPEAT_THRESHOLD_SET]
                      [-rl] [-mqt MAPQ_THRESHOLD] [-d]
                      [-r REPEATS [REPEATS ...]] [-con] [-gc1 LOWERGC]
                      [-gc2 UPPERGC] [-nf]
                      [-rc TVRS_FOR_CONTEXT [TVRS_FOR_CONTEXT ...]]
                      [-bp BP_CONTEXT] [-pl] [-pff {pdf,png,svg,all}] [-p1]
                      [-p2] [-p3] [-p4] [-p5] [-p6] [-p7] [-p8] [-prc]

Estimation of telomere content from WGS data of a tumor and/or a control
sample.

optional arguments:
  -h, --help            show this help message and exit
  -ibt TUMOR_BAM, --inputBamTumor TUMOR_BAM
                        Path to the indexed input BAM file of the tumor
                        sample.
  -ibc CONTROL_BAM, --inputBamControl CONTROL_BAM
                        Path to the indexed input BAM file of the control
                        sample.
  -o OUTPUT_DIR, --outPath OUTPUT_DIR
                        Path to the output directory into which all results
                        are written.
  -p PID, --pid PID     Sample name used in output files and diagrams
                        (required).
  -b BANDING_FILE, --bandingFile BANDING_FILE
                        Path to a tab-separated file with information on
                        chromosome banding. The first four columns of the
                        table have to contain the chromosome name, the start
                        and end position and the band name. The table should
                        not have a header. If no banding file is specified,
                        the banding information of hg19 will be used.
  -rt REPEAT_THRESHOLD_SET, --repeatThreshold REPEAT_THRESHOLD_SET
                        The number of repeats needed for a read to be
                        classified as telomeric. If no repeat threshold is
                        defined, TelomereHunter will calculate the
                        repeat_threshold depending on the read length with the
                        following formula: repeat_threshold =
                        floor(read_length * 6/100)
  -rl, --perReadLength  Repeat threshold is set per 100 bp read length. The
                        used repeat threshold will be: floor(read_length *
                        repeat_threshold/100) E.g. Setting -rt 8 -rl means
                        that 8 telomere repeats are required per 100 bp read
                        length. If the read length is 50 bp, the threshold is
                        set to 4.
  -mqt MAPQ_THRESHOLD, --mappingQualityThreshold MAPQ_THRESHOLD
                        The mapping quality needed for a read to be considered
                        as mapped (default = 8).
  -d, --removeDuplicates
                        Reads marked as duplicates in the input bam file(s)
                        are removed in the filtering step.
  -r REPEATS [REPEATS ...], --repeats REPEATS [REPEATS ...]
                        List of telomere repeat types to search for. Reverse
                        complements are automatically generated and do not
                        need to be specified! By default, TelomereHunter
                        searches for t-, g-, c- and j-type repeats (TTAGGG
                        TGAGGG TCAGGG TTGGGG).
  -con, --consecutive   Search for consecutive repeats.
  -gc1 LOWERGC, --lowerGC LOWERGC
                        Lower limit used for GC correction of telomere
                        content. The value must be an integer between 0 and
                        100 (default = 48).
  -gc2 UPPERGC, --upperGC UPPERGC
                        Upper limit used for GC correction of telomere
                        content. The value must be an integer between 0 and
                        100 (default = 52).
  -nf, --noFiltering    If the filtering step of TelomereHunter has already
                        been run previously, skip this step.
  -rc TVRS_FOR_CONTEXT [TVRS_FOR_CONTEXT ...], --repeatsContext TVRS_FOR_CONTEXT [TVRS_FOR_CONTEXT ...]
                        List of telomere variant repeats for which to analyze
                        the sequence context. Reverse complements are
                        automatically generated and do not need to be
                        specified! Counts for these telomere variant repeats
                        (arbitrary and singleton context) will be added to the
                        summary table. Default repeats: TCAGGG TGAGGG TTGGGG
                        TTCGGG TTTGGG ATAGGG CATGGG CTAGGG GTAGGG TAAGGG).
  -bp BP_CONTEXT, --bpContext BP_CONTEXT
                        Number of base pairs on either side of the telomere
                        variant repeat to investigate. Please use a number
                        that is divisible by 6.
  -pl, --parallel       The filtering, sorting and estimating steps of the
                        tumor and control sample are run in parallel. This
                        will speed up the computation time of TelomereHunter.
  -pff {pdf,png,svg,all}, --plotFileFormat {pdf,png,svg,all}
                        File format of output diagrams. Choose from pdf
                        (default), png, svg or all (pdf, png and svg).
  -p1, --plotChr        Make diagrams with telomeric reads mapping to each
                        chromosome. If none of the options p1/p2/p3/p4/p5/p6
                        are chosen, all diagrams will be created.
  -p2, --plotFractions  Make a diagram with telomeric reads in each fraction
                        (intrachromosomal, subtelomeric, junction spanning,
                        intratelomeric). If none of the options
                        p1/p2/p3/p4/p5/p6 are chosen, all diagrams will be
                        created.
  -p3, --plotTelContent
                        Make a diagram with the gc corrected telomere content
                        in the analyzed samples. If none of the options
                        p1/p2/p3/p4/p5/p6 are chosen, all diagrams will be
                        created.
  -p4, --plotGC         Make a diagram with GC content distributions in all
                        reads and in intratelomeric reads. If none of the
                        options p1/p2/p3/p4/p5/p6 are chosen, all diagrams
                        will be created.
  -p5, --plotRepeatFreq
                        Make histograms of the repeat frequencies per
                        intratelomeric read. If none of the options
                        p1/p2/p3/p4/p5/p6 are chosen, all diagrams will be
                        created.
  -p6, --plotTVR        Make plots for telomere variant repeats.
  -p7, --plotSingleton  Make plots for singleton telomere variant repeats.
  -p8, --plotNone       Do not make any diagrams. If none of the options
                        p1/p2/p3/p4/p5/p6/p7/p8 are chosen, all diagrams will
                        be created.
  -prc, --plotRevCompl  Distinguish between forward and reverse complement
                        telomere repeats in diagrams.

Contact Lina Sieverling (l.sieverling@dkfz-heidelberg.de) for questions and
support.

[user@cn3200 ~]$ **mkdir output\_dir**

[user@cn3200 ~]$ **telomerehunter -o output\_dir -p test\_sample -ibt $TH\_DATA/test\_sample\_control.bam -ibc $TH\_DATA/test\_sample\_tumor.bam -p test\_sample -b $TH\_DATA/hg38\_cytoband\_telomerehunter.txt -pl 2> test\_sample.err.txt** 


        TelomereHunter  Copyright (C) 2015  Lina Sieverling, Philip Ginsbach, Lars Feuerbach
        This program comes with ABSOLUTELY NO WARRANTY.
        This is free software, and you are welcome to redistribute it
        under certain conditions. For details see the GNU General Public License
        in the license copy received with TelomereHunter or .


TelomereHunter 1.1.0



Repeat threshold was not set by user. Setting it to 6 reads per 100 bp read length.
Calculating repeat threshold for tumor sample: Read length is 151. Repeat threshold is set to 9.
Calculating repeat threshold for control sample: Read length is 151. Repeat threshold is set to 9.

------ Tumor Sample: started filtering telomere reads ------
------ Control Sample: started filtering telomere reads ------
...

```

End the interactive session:

```

[user@cn3200 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





