

document.querySelector('title').textContent = 'OncodriveFML: identifying coding and non-coding regions with cancer driver mutations';
**OncodriveFML: identifying coding and non-coding regions with cancer driver mutations**


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



OncodriveFML is a method designed to analyze the pattern of somatic mutations 
across tumors in both coding and non-coding genomic regions 
to identify signals of positive selection, and therefore, their involvement in tumorigenesis.
It can be used to identify protein-coding genes, promoters, untranslated regions,
intronic splice regions, and lncRNAs-containing driver mutations.



### References:


* Loris Mularoni, Radhakrishnan Sabarinathan, Jordi Deu-Pons, Abel Gonzalez-Perez and Núria López-Biga,   

*OncodriveFML: a general framework to identify coding and non-coding regions with cancer driver mutations*   

[Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0994-0), 2016 **17**-128. DOI 10.1186/s13059-016-0994-0


Documentation
* [OncodriveFML Bitbucket repository](https://bitbucket.org/bbglab/oncodrivefml/src/master/)
* [OncodriveFML Documentation](https://oncodrivefml.readthedocs.io/en/latest/oncodriveFML.html)


Important Notes
* Module Name: oncodriveFML (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **ODFML\_HOME**  installation directory
	+ **ODFML\_BIN**       executable directory
	+ **ODFML\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
[user@cn3335 ~]$ **module load oncodriveFML**
[+] Loading singularity  3.10.5  on cn4338
[+] Loading oncodriveFML  2.2.0
[user@cn3335 ~]$ **oncodrivefml -h**
Usage: oncodrivefml [OPTIONS]

  Run OncodriveFML on the genomic regions in ELEMENTS FILE using the mutations
  in MUTATIONS FILE.

Options:
  -i, --input MUTATIONS_FILE      Variants file  [required]
  -e, --elements ELEMENTS_FILE    Genomic elements to analyse  [required]
  -t, --type [coding|noncoding]   Deprecated option
  -s, --sequencing [wgs|wes|targeted]
                                  Type of sequencing: whole genome, whole
                                  exome or targeted.
  -o, --output OUTPUT_FOLDER      Output folder. Default to regions file name
                                  without extensions.
  -c, --configuration CONFIG_FILE
                                  Configuration file. Default to
                                  'oncodrivefml_v2.conf' in the current folder
                                  if exists or to
                                  ~/.config/bbglab/oncodrivefml_v2.conf if
                                  not.
  --samples-blacklist SAMPLES_BLACKLIST
                                  Remove these samples when loading the input
                                  file.
  --signature SIGNATURE           File with the signatures to use
  --signature-correction [wg|wx]  Correct the computed signutares by genomic
                                  or exomic signtures. Only valid for human
                                  genomes
  --no-indels                     Discard indels in your analysis
  --cores INTEGER                 Cores to use. Default: all
  --seed INTEGER RANGE            Set up an initial random seed to have
                                  reproducible results  [0<x<=4294967295]
  --generate-pickle               Deprecated flag. Do not use.
  --force                         Overwrite results if exists
  --debug                         Show more progress details
  --version                       Show the version and exit.
  -h, --help                      Show this message and exit.

```

Copy sample data to the current folder:

```

[user@cn3335 ~]$ **cp $ODFML\_DATA/\* .**

```

Run oncodriveFML on the sample data:

```

[user@cn3335 ~]$ **oncodrivefml -i paad.txt.gz -e cds.tsv.gz --sequencing wes**
...

```





