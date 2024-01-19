

document.querySelector('title').textContent = 'funannotate: a pipeline for genome annotation ';
**funannotate: a pipeline for genome annotation** 


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



Funannotate is a genome prediction, annotation, and comparison software package. 
It was originally written to annotate fungal genomes (small eukaryotes ~ 30 Mb genomes), 
but has evolved over time to accomodate larger genomes. 



### References:


* Palmer J   

 [*Funannotate: Fungal genome annotation scripts (2017).*](https://github.com/nextgenusfs/funannotate)
* Wan-Chen Li and Ting-Fang Wang   

*PacBio Long-Read Sequencing, Assembly, and Funannotate Reannotation of the Complete Genome of Trichoderma reesei QM6a.*   

[Chapter 21 in: Trichoderma reesei, Methods in Molecular Biology, vol. 2234, pp.311-329](https://link.springer.com/protocol/10.1007/978-1-0716-1048-0_21)


Documentation
* [Funannotate documentation](https://funannotate.readthedocs.io/en/latest/)
* [Funannotate GitHub page](https://github.com/nextgenusfs/funannotate)


Important Notes
* Module Name: Funannotate (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **FA\_HOME**  installation directory
	+ **FA\_BIN**    executables directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
  


```

[user@biowulf]$ **sinteractive --gres=lscratch:10 --mem=10g -c4**
[user@cn0861 ~]$ **module load funannotate** 
[+] Loading funannotate  1.8.9  ...      
[+] Loading signalp 4.1  ...
[+] Loading singularity  3.8.0  ...      

```

The funannotate data processing pipeline involves several core modules:

```

[user@biowulf]$ **funannotate -h** 

Usage:       funannotate <command> <arguments>
version:     1.8.9

Description: Funannotate is a genome prediction, annotation, and comparison pipeline.

Commands:
  clean       Find/remove small repetitive contigs
  sort        Sort by size and rename contig headers
  mask        Repeatmask genome assembly

  train       RNA-seq mediated training of Augustus/GeneMark
  predict     Run gene prediction pipeline
  fix         Fix annotation errors (generate new GenBank file)
  update      RNA-seq/PASA mediated gene model refinement
  remote      Partial functional annotation using remote servers
  iprscan     InterProScan5 search (Docker or local)
  annotate    Assign functional annotation to gene predictions
  compare     Compare funannotated genomes

  util        Format conversion and misc utilities
  setup       Setup/Install databases
  test        Download/Run funannotate installation tests
  check       Check Python, Perl, and External dependencies [--show-versions]
  species     list pre-trained Augustus species
  database    Manage databases
  outgroups   Manage outgroups for funannotate compare

```

These modules are supposed to be used in the following order:   


clean → sort → mask → train → predict → update → annotate → compare   


```


[user@biowulf]$ **funannotate setup -h**

Usage:       funannotate clean <arguments>
version:     1.8.9

Description: The script sorts contigs by size, starting with shortest contigs it uses minimap2
             to find contigs duplicated elsewhere, and then removes duplicated contigs.

Arguments:
  -i, --input    Multi-fasta genome file (Required)
  -o, --out      Cleaned multi-fasta output file (Required)
  -p, --pident   Percent identity of overlap. Default = 95
  -c, --cov      Percent coverage of overlap. Default = 95
  -m, --minlen   Minimum length of contig to keep. Default = 500
  --exhaustive   Test every contig. Default is to stop at N50 value.

[user@biowulf]$ **funannotate setup -i pfam -d ./databases --force** 
-------------------------------------------------------
[Aug 15 07:50 PM]: OS: Debian GNU/Linux 10, 56 cores, ~ 264 GB RAM. Python: 3.8.13
[Aug 15 07:50 PM]: Running 1.8.13
[Aug 15 07:50 PM]: Database location: ./databases
[Aug 15 07:50 PM]: Retrieving download links from GitHub Repo
[Aug 15 07:50 PM]: Parsing Augustus pre-trained species and porting to funannotate
[Aug 15 07:50 PM]: Downloading UniProtKB/SwissProt database
[Aug 15 07:50 PM]: Downloading: http://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz Bytes: 91192930
[Aug 15 07:50 PM]: Downloading: http://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt Bytes: 151
[Aug 15 07:50 PM]: Building diamond database
[Aug 15 07:50 PM]: UniProtKB Database: version=2022_03 date=2022-08-03 records=568,002
[Aug 15 07:50 PM]: Funannoate setup complete. Add this to ~/.bash_profile or ~/.bash_aliases:

        export FUNANNOTATE_DB=/gpfs/gsfs7/users/user/funannotate/databases

[user@biowulf]$ **funannotate sort -h**

Usage:       funannotate sort <arguments>
version:     1.8.9

Description: Script will download/format necessary databases for funannotate.

Options:
  -i, --install    Download format databases. Default: all
                     [merops,uniprot,dbCAN,pfam,repeats,go,
                      mibig,interpro,busco_outgroups,gene2product]
  -b, --busco_db   Busco Databases to install. Default: dikarya [all,fungi,aves,etc]
  -d, --database   Path to funannotate database
  -u, --update     Check remote md5 and update if newer version found
  -f, --force      Force overwriting database
  -w, --wget       Use wget to download instead of python requests
  -l, --local      Use local resource JSON file instead of current on github
[user@biowulf]$ **funannotate -i pfam** 
[Aug 15 01:23 PM]: OS: CentOS Linux 7, 56 cores, ~ 264 GB RAM. Python: 3.7.11
[Aug 15 01:23 PM]: Running 1.8.9
[Aug 15 01:23 PM]: Database location: /fdb/funannotate/db
[Aug 15 01:23 PM]: Retrieving download links from GitHub Repo
[Aug 15 01:23 PM]: Parsing Augustus pre-trained species and porting to funannotate
[Aug 15 01:23 PM]: Pfam Database: version=34.0 date=2021-03 records=19,179
[user@biowulf]$ **funannotate mask -h** 

Usage:       funannotate mask <arguments>
version:     1.8.9

Description: This script is a wrapper for repeat masking. Default is to run very simple
             repeat masking with tantan. The script can also run RepeatMasker and/or
             RepeatModeler. It will generate a softmasked genome. Tantan is probably not
             sufficient for soft-masking an assembly, but with RepBase no longer being
             available RepeatMasker/Modeler may not be functional for many users.

Arguments:
  -i, --input                  Multi-FASTA genome file. (Required)
  -o, --out                    Output softmasked FASTA file. (Required)

Optional:
  -m, --method                 Method to use. Default: tantan [repeatmasker, repeatmodeler]
  -s, --repeatmasker_species   Species to use for RepeatMasker
  -l, --repeatmodeler_lib      Custom repeat database (FASTA format)
  --cpus                       Number of cpus to use. Default: 2
  --debug                      Keep intermediate files

```

etc.


```

[user@cn0861 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226

```





