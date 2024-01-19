

document.querySelector('title').textContent = 'VADR: Viral Annotation DefineR.';
VADR: Viral Annotation DefineR.


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



VADR is a suite of tools for classifying and analyzing sequences 
homologous to a set of reference models of viral genomes or gene families.
It has been mainly tested for analysis of Norovirus, Dengue, and SARS-CoV-2 virus sequences 
in preparation for submission to the GenBank database.



### References:


* Alejandro A Schäffer, Eneida L Hatcher, Linda Yankie, Lara Shonkwiler, J Rodney Brister, Ilene Karsch-Mizrachi, Eric P Nawrocki   

 *VADR: validation and annotation of virus sequence submissions to GenBank`.*  

BMC Bioinformatics (2020) 21:211. https://doi.org/10.1186/s12859-020-3537-3


Documentation
* [VADR Github page](https://github.com/ncbi/vadr)


Important Notes
* Module Name: VADR (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **VADR\_HOME**  VADR installation directory
	+ **VADR\_BIN**       VADR executable directory
	+ **VADR\_SRC**     a folder containing the source code
	+ **VADR\_DATA**    a folder containing sample data
	+ **VADR\_HMM**     a folder containing sample Hidden Markov Models



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive** 
[user@cn4471 ~]$ **module load vadr** 
[+] Loading hmmer  3.3.2  on cn0847
[+] Loading perl 5.24.3 on cn0847
[+] Loading vadr / 1.1.3  ...

```

Copy sample data to your current folder:

```

[user@cn4471 ~]$ **cp -r $VADR\_DATA/\* .**

```

Preprocess the data with samtools:

```

[user@cn4471 ~]$ **v-build.pl -h**
...
# v-build.pl :: build homology model of a single sequence for feature annotation
# VADR 1.1.3 (Feb 2021)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:    Mon Mar 29 12:40:22 2021
#
Usage: v-build.pl [-options] <accession> <path to output directory to create>

basic options:
  -f             : force; if dir <output directory> exists, overwrite it
  -v             : be verbose; output commands to stdout as they're run
  --stk <s>      : read single sequence stockholm 'alignment' from <s>
  --infa <s>     : read single sequence fasta file from <s>, don't fetch it
  --inft <s>     : read feature table file from <s>, don't fetch it
  --ftfetch1     : fetch feature table with efetch -format ft
  --ftfetch2     : fetch feature table with efetch -format gbc | xml2tbl
  --gb           : parse a genbank file, not a feature table file
  --ingb <s>     : read genbank file from <s>, don't fetch it
  --addminfo <s> : add feature info from model info file <s>
  --forcelong    : allow long models > 25Kb in length
  --keep         : do not remove intermediate files, keep them all on disk

options for controlling what feature types are stored in model info file
[default set is: CDS,gene,mat_peptide]:
  --fall      : store info for all feature types (except those in --fskip)
  --fadd <s>  : also store feature types in comma separated string <s>
  --fskip <s> : do not store info for feature types in comma separated string <s>

options for controlling what qualifiers are stored in model info file
[default set is:product,gene,exception]:
  --qall        : store info for all qualifiers (except those in --qskip)
  --qadd <s>    : also store info for qualifiers in comma separated string <s>
  --qftradd <s> : --qadd <s2> only applies for feature types in comma separated string <s>
  --qskip <s>   : do not store info for qualifiers in comma separated string <s>
  --noaddgene   : do not add gene qualifiers from gene features to overlapping features

options for including additional model attributes:
  --group <s>    : specify model group is <s>
  --subgroup <s> : specify model subgroup is <s>

options for controlling CDS translation step:
  --ttbl <n> : use NCBI translation table <n> to translate CDS [1]

options for controlling cmbuild step:
  --cmn <n>       : set number of seqs for glocal fwd HMM calibration to <n>
  --cmp7ml        : set CM's filter p7 HMM as the ML p7 HMM
  --cmere <x>     : set CM relative entropy target to <x>
  --cmeset <x>    : set CM eff seq # for CM to <x>
  --cmemaxseq <x> : set CM maximum alowed eff seq # for CM to <x>
  --cminfile <s>  : read cmbuild options from file <s>

options for skipping stages:
  --skipbuild : skip the cmbuild and blastn db creation steps
  --onlyurl   : output genbank file url for accession and exit

optional output files:
  --ftrinfo : create file with internal feature information
  --sgminfo : create file with internal segment information

other expert options:
  --execname <s> : define executable name of this script as <s>

```


```

[user@cn4471 ~]$ **v-build.pl NC\_039897 NC\_039897**
# v-build.pl :: build homology model of a single sequence for feature annotation
# VADR 1.1.3 (Feb 2021)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:              Tue Mar 30 14:51:04 2021
# $VADRBLASTDIR:     /usr/local/apps/VADR/1.1.3/vadr-vadr-1.1.3/ncbi-blast/bin
# $VADREASELDIR:     /usr/local/apps/VADR/1.1.3/vadr-vadr-1.1.3/Bio-Easel/src/easel/miniapps
# $VADRINFERNALDIR:  /usr/local/apps/VADR/1.1.3/vadr-vadr-1.1.3/infernal/binaries
# $VADRSCRIPTSDIR:   /usr/local/apps/VADR/1.1.3/vadr-vadr-1.1.3/vadr
#
# accession/model name:  NC_039897
# output directory:      NC_039897
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fetching FASTA file                                          ... done. [    3.3 seconds]
# Parsing FASTA file                                           ... done. [    0.0 seconds]
# Fetching feature table file                                  ... done. [    3.9 seconds]
# Parsing feature table file                                   ... done. [    0.0 seconds]
# Fetching and parsing protein feature table file(s)           ... done. [   11.8 seconds]
# Pruning data read from GenBank                               ... done. [    0.0 seconds]
# Reformatting FASTA file to Stockholm file                    ... done. [    0.0 seconds]
# Finalizing feature information                               ... done. [    0.0 seconds]
# Translating CDS                                              ... done. [    0.0 seconds]
# Building BLAST protein database                              ... done. [    0.2 seconds]
# Building HMMER protein database                              ... done. [    2.6 seconds]
# Building CM (should take roughly 10-30 minutes)              ... done. [  749.1 seconds]
# Pressing CM file                                             ... done. [    0.2 seconds]
# Building BLAST nucleotide database of CM consensus           ... done. [    0.4 seconds]
# Creating model info file                                     ... done. [    0.0 seconds]
#
# Output printed to screen saved in:                                           NC_039897.vadr.log
# List of executed commands saved in:                                          NC_039897.vadr.cmd
# List and description of all output files saved in:                           NC_039897.vadr.filelist
# fasta file for NC_039897 saved in:                                           NC_039897.vadr.fa
# feature table format file for NC_039897 saved in:                            NC_039897.vadr.tbl
# feature table format file for YP_009538340.1 saved in:                       NC_039897.vadr.YP_009538340.1.tbl
# feature table format file for YP_009538341.1 saved in:                       NC_039897.vadr.YP_009538341.1.tbl
# feature table format file for YP_009538342.1 saved in:                       NC_039897.vadr.YP_009538342.1.tbl
# Stockholm alignment file for NC_039897 saved in:                             NC_039897.vadr.stk
# fasta sequence file for CDS from NC_039897 saved in:                         NC_039897.vadr.cds.fa
# fasta sequence file for translated CDS from NC_039897 saved in:              NC_039897.vadr.protein.fa
# BLAST db .phr file for NC_039897 saved in:                                   NC_039897.vadr.protein.fa.phr
# BLAST db .pin file for NC_039897 saved in:                                   NC_039897.vadr.protein.fa.pin
# BLAST db .psq file for NC_039897 saved in:                                   NC_039897.vadr.protein.fa.psq
# BLAST db .pdb file for NC_039897 saved in:                                   NC_039897.vadr.protein.fa.pdb
# BLAST db .pot file for NC_039897 saved in:                                   NC_039897.vadr.protein.fa.pot
# BLAST db .ptf file for NC_039897 saved in:                                   NC_039897.vadr.protein.fa.ptf
# BLAST db .pto file for NC_039897 saved in:                                   NC_039897.vadr.protein.fa.pto
# HMMER model db file for NC_039897 saved in:                                  NC_039897.vadr.protein.hmm
# hmmbuild build output (concatenated) saved in:                               NC_039897.vadr.protein.hmmbuild
# binary HMM and p7 HMM filter file saved in:                                  NC_039897.vadr.protein.hmm.h3m
# SSI index for binary HMM file saved in:                                      NC_039897.vadr.protein.hmm.h3i
# optimized p7 HMM filters (MSV part) saved in:                                NC_039897.vadr.protein.hmm.h3f
# optimized p7 HMM filters (remainder) saved in:                               NC_039897.vadr.protein.hmm.h3p
# hmmpress output file saved in:                                               NC_039897.vadr.hmmpress
# CM file saved in:                                                            NC_039897.vadr.cm
# cmbuild output file saved in:                                                NC_039897.vadr.cmbuild
# binary CM and p7 HMM filter file saved in:                                   NC_039897.vadr.cm.i1m
# SSI index for binary CM file saved in:                                       NC_039897.vadr.cm.i1i
# optimized p7 HMM filters (MSV part) saved in:                                NC_039897.vadr.cm.i1f
# optimized p7 HMM filters (remainder) saved in:                               NC_039897.vadr.cm.i1p
# cmpress output file saved in:                                                NC_039897.vadr.cmpress
# fasta sequence file with cmemit consensus sequence for NC_039897 saved in:   NC_039897.vadr.nt.fa
# BLAST db .nhr file for NC_039897 saved in:                                   NC_039897.vadr.nt.fa.nhr
# BLAST db .nin file for NC_039897 saved in:                                   NC_039897.vadr.nt.fa.nin
# BLAST db .nsq file for NC_039897 saved in:                                   NC_039897.vadr.nt.fa.nsq
# BLAST db .ndb file for NC_039897 saved in:                                   NC_039897.vadr.nt.fa.ndb
# BLAST db .not file for NC_039897 saved in:                                   NC_039897.vadr.nt.fa.not
# BLAST db .ntf file for NC_039897 saved in:                                   NC_039897.vadr.nt.fa.ntf
# BLAST db .nto file for NC_039897 saved in:                                   NC_039897.vadr.nt.fa.nto
# VADR 'model info' format file for NC_039897 saved in:                        NC_039897.vadr.minfo
#
# All output files created in directory ./NC_039897/
#
# Elapsed time:  00:12:51.55
#                hh:mm:ss
#
[ok]

```


```

[user@cn4471 ~]$ **v-annotate.pl -h** 
# v-annotate.pl :: classify and annotate sequences using a CM library
# VADR 1.1.3 (Feb 2021)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:    Tue Mar 30 15:17:51 2021
#
Usage: v-annotate.pl [-options] <fasta file to annotate> <output directory to create>

basic options:
  -f             : force; if output dir exists, overwrite it
  -v             : be verbose; output commands to stdout as they're run
  --atgonly      : only consider ATG a valid start codon
  --minpvlen <n> : min CDS/mat_peptide/gene length for feature table output and protein validation is <n> [30]
  --keep         : do not remove intermediate files, keep them all on disk

options for specifying classification:
  --group <s>    : set expected classification of all seqs to group <s>
  --subgroup <s> : set expected classification of all seqs to subgroup <s>

options for controlling severity of alerts:
  --alt_list        : output summary of all alerts and exit
  --alt_pass <s>    : specify that alert codes in comma-separated <s> do not cause FAILure
  --alt_fail <s>    : specify that alert codes in comma-separated <s> do cause FAILure
  --alt_mnf_yes <s> : alert codes in <s> for 'misc_not_failure' features cause misc_feature-ization, not failure
  --alt_mnf_no <s>  : alert codes in <s> for 'misc_not_failure' features cause failure, not misc-feature-ization
  --ignore_mnf      : ignore non-zero 'misc_not_feature' values in .minfo file, set to 0 for all features/models

options related to model files:
  -m <s>      : use CM file <s> instead of default
  -a <s>      : use HMM file <s> instead of default
  -i <s>      : use model info file <s> instead of default
  -n <s>      : use blastn db file <s> instead of default
  -x <s>      : blastx dbs are in dir <s>, instead of default
  --mkey <s>  : .cm, .minfo, blastn .fa files in $VADRMODELDIR start with key <s>, not 'vadr'
  --mdir <s>  : model files are in directory <s>, not in $VADRMODELDIR
  --mlist <s> : only use models listed in file <s>

options for controlling output feature table:
  --nomisc        : in feature table for failed seqs, never change feature type to misc_feature
  --notrim        : in feature table, don't trim coords due to Ns (for any feature types)
  --noftrtrim <s> : in feature table, don't trim coords due to Ns for feature types in comma-delmited <s>
  --noprotid      : in feature table, don't add protein_id for CDS and mat_peptides
  --forceprotid   : in feature table, force protein_id value to be sequence name, then idx

options for controlling thresholds related to alerts:
  --lowsc <x>       : lowscore/LOW_SCORE bits per nucleotide threshold is <x> [0.3]
  --indefclass <x>  : indfcls/INDEFINITE_CLASSIFICATION bits per nucleotide diff threshold is <x> [0.03]
  --incspec <x>     : inc{group,subgrp}/INCORRECT_{GROUP,SUBGROUP} bits/nt threshold is <x> [0.2]
  --lowcov <x>      : lowcovrg/LOW_COVERAGE fractional coverage threshold is <x> [0.9]
  --dupregolp <n>   : dupregin/DUPLICATE_REGIONS minimum model overlap is <n> [20]
  --dupregsc <x>    : dupregin/DUPLICATE_REGIONS minimum bit score is <x> [10]
  --indefstr <x>    : indfstrn/INDEFINITE_STRAND minimum weaker strand bit score is <x> [25]
  --lowsim5term <n> : lowsim5{s,f}/LOW_{FEATURE_}SIMILARITY_START minimum length is <n> [15]
  --lowsim3term <n> : lowsim3{s,f}/LOW_{FEATURE_}SIMILARITY_END minimum length is <n> [15]
  --lowsimint <n>   : lowsimi{s,f}/LOW_{FEATURE_}SIMILARITY (internal) minimum length is <n> [1]
  --biasfract <x>   : biasdseq/BIASED_SEQUENCE fractional threshold is <x> [0.25]
  --indefann <x>    : indf{5,3}loc/'INDEFINITE_ANNOTATION_{START,END} non-mat_peptide min allowed post probability is <x> [0.8]
  --indefann_mp <x> : indf{5,3}loc/'INDEFINITE_ANNOTATION_{START,END} mat_peptide min allowed post probability is <x> [0.6]
  --fstminnt <n>    : fst{hi,lo}cnf/POSSIBLE_FRAMESHIFT_{HIGH,LOW}_CONF max allowed frame disagreement nt length w/o alert is <n> [6]
  --fsthighthr <x>  : fsthicnf/POSSIBLE_FRAMESHIFT_HIGH_CONF minimum average probability for alert is <x> [0.8]
  --fstlowthr <x>   : fstlocnf/POSSIBLE_FRAMESHIFT_LOW_CONF minimum average probability for alert is <x> [0.3]
  --xalntol <n>     : indf{5,3}{st,lg}/INDEFINITE_ANNOTATION_{START,END} max allowed nt diff blastx start/end is <n> [5]
  --xmaxins <n>     : insertnp/INSERTION_OF_NT max allowed nucleotide insertion length in blastx validation is <n> [27]
  --xmaxdel <n>     : deletinp/DELETION_OF_NT max allowed nucleotide deletion length in blastx validation is <n> [27]
  --nmaxins <n>     : insertnn/INSERTION_OF_NT max allowed nucleotide (nt) insertion length in CDS nt alignment is <n> [27]
  --nmaxdel <n>     : deletinn/DELETION_OF_NT max allowed nucleotide (nt) deletion length in CDS nt alignment is <n> [27]
  --xlonescore <n>  : indfantp/INDEFINITE_ANNOTATION min score for a blastx hit not supported by CM analysis is <n> [80]
  --hlonescore <n>  : indfantp/INDEFINITE_ANNOTATION min score for a hmmer hit not supported by CM analysis is <n> [10]

options for controlling cmalign alignment stage:
  --mxsize <n> : set max allowed memory for cmalign to <n> Mb [16000]
  --tau <x>    : set the initial tau value for cmalign to <x> [0.001]
  --nofixedtau : do not fix the tau value when running cmalign, allow it to decrease if nec
  --nosub      : use alternative alignment strategy for truncated sequences
  --noglocal   : do not run cmalign in glocal mode (run in local mode)

options for controlling blastx protein validation stage:
  --xmatrix <s> : use the matrix <s> with blastx (e.g. BLOSUM45)
  --xdrop <n>   : set the xdrop value for blastx to <n> [25]
  --xnumali <n> : number of alignments to keep in blastx output and consider if --xlongest is <n> [20]
  --xlongest    : keep the longest blastx hit, not the highest scoring one

options for using hmmer instead of blastx for protein validation:
  --hmmer        : use hmmer for protein validation, not blastx
  --h_max        : use --max option with hmmsearch
  --h_minbit <x> : set minimum hmmsearch bit score threshold to <x> [-10]

options related to blastn-derived seeded alignment acceleration:
  -s               : use the max length ungapped region from blastn to seed the alignment
  --s_blastnws <n> : for -s, set blastn -word_size <n> to <n> [7]
  --s_blastnsc <x> : for -s, set blastn minimum HSP score to consider to <x> [50]
  --s_overhang <n> : for -s, set length of nt overhang for subseqs to align to <n> [100]

options related to replacing Ns with expected nucleotides:
  -r               : replace stretches of Ns with expected nts, where possible
  --r_minlen <n>   : minimum length subsequence to replace Ns in is <n> [5]
  --r_minfract <x> : minimum fraction of Ns in subseq to trigger replacement is <x> [0.5]
  --r_fetchr       : fetch features for output fastas from seqs w/Ns replaced, not originals
  --r_cdsmpr       : detect CDS and MP alerts in sequences w/Ns replaced, not originals
  --r_pvorig       : use original sequences for protein validation, not replaced seqs
  --r_prof         : use slower profile methods, not blastn, to identify Ns to replace

options related to parallelization on compute farm:
  -p             : parallelize cmsearch/cmalign on a compute farm
  -q <s>         : use qsub info file <s> instead of default
  --nkb <n>      : number of KB of sequence for each farm job is <n> [10]
  --wait <n>     : allow <n> wall-clock minutes for jobs on farm to finish, including queueing time [500]
  --errcheck     : consider any farm stderr output as indicating a job failure
  --maxnjobs <n> : set max number of jobs to submit to compute farm to <n> [2500]

options for skipping stages:
  --skip_align : skip the cmalign step, use results from an earlier run of the script
  --skip_pv    : do not perform blastx-based protein validation

optional output files:
  --out_stk     : output per-model full length stockholm alignments (.stk)
  --out_afa     : output per-model full length fasta alignments (.afa)
  --out_rpstk   : with -r, output stockholm alignments of seqs with Ns replaced
  --out_rpafa   : with -r, output fasta alignments of seqs with Ns replaced
  --out_nofs    : do not output frameshift stockholm alignment files
  --out_nofasta : do not output fasta files of features, or passing/failing seqs
  --out_debug   : dump voluminous info from various data structures to output files

other expert options:
  --execname <s> : define executable name of this script as <s>
  --alicheck     : for debugging, check aligned sequence vs input sequence for identity
  --noseqnamemax : do not enforce a maximum length of 50 for sequence names (GenBank max)
  --minbit <x>   : set minimum cmsearch bit score threshold to <x> [-10]
  --origfa       : do not copy fasta file prior to analysis, use original
  --msub <s>     : read model substitution file from <s>
  --xsub <s>     : read blastx db substitution file from <s>

```

End the interactive session:

```

[user@cn4471 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





