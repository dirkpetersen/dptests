

document.querySelector('title').textContent = 'ANNOgesic: Accurate RNA-Seq-based annotation of bacterial and archaeal genomes';
**ANNOgesic: Accurate RNA-Seq-based annotation of bacterial and archaeal genomes**


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



Processing and integrating RNA-Seq data in order to
generate high-resolution annotations is challenging, time consuming and requires
numerous different steps. ANNOgesic is a powerful and modular
pipeline that provides the required analyses and simplifies
RNA-Seq-based bacterial and archaeal genome annotation. It predicts and
annotates numerous features, including small non-coding RNAs, with high
precision.



### References:


* Sung-Huan Yu, Jörg Vogel and Konrad U.Förstner, 2018.  

*ANNOgesic: a Swiss army knife for the RNA-seq based
annotation of bacterial/archaeal genomes.*   

[GigaScience, DOI:10.1093/gigascience/giy096, PMID:30169674.](https://academic.oup.com/gigascience/article/7/9/giy096/5087959)* Sung-Huan Yu, Jörg Vogel and Konrad Ulrich Förstner,   

*ANNOgesic: A Pipeline To Translate Bacterial/Archaeal RNA-Seq Data Into High-Resolution Genome Annotations.*   
[bioRxiv preprint](https://www.biorxiv.org/content/early/2017/05/29/143081) first posted online May. 27, 2017; doi: http://dx.doi.org/10.1101/143081i.


Documentation
* [ANNOgesic documentation](https://annogesic.readthedocs.io/en/latest/)
* [ANNOgesic project page](https://pypi.org/project/ANNOgesic/)
* [ANNOgesic Github page](https://github.com/Sung-Huan/ANNOgesic)


Important Notes
* Module Name: ANNOgesic (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Implemented as a Singularity container
* Unusual environment variables set
	+ **AG\_HOME**  ANNOgesic installation directory
	+ **AG\_BIN**       ANNOgesic executable directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=4g**
[user@cn3316 ~]$ **module load ANNOgesic**
[+] Loading singularity  on cn3316 
[+] Loading ANNOgesic 1.0.2  ...

```

At this point, user has two options:
  
 1) typing the command

```

[user@cn3316 user]$ **ag** 

```

(without arguments) will bring the user into the singularity container shell environment

```

Singularity ANNOgesic.sqsh:~>

```

from which one can run any script or command accessible within the container on any data accessible from inside the container. For example, the following commands will run built-in tests:

```

Singularity ANNOgesic.sqsh:~> **python /ANNOgesic/tests/test\_operon.py** 
Detecting operons of test
Warning: No proper file - test.gff
..
----------------------------------------------------------------------
Ran 2 tests in 0.009s
  

Singularity ANNOgesic.sqsh:~> **python3 /ANNOgesic/tests/test\_plot\_PPI.py** 
.......Plotting nusB
..
----------------------------------------------------------------------
Ran 9 tests in 1.567s

OK

```

To exit from the container shell environment, type:

```
Singularity ANNOgesic.sqsh:~> **exit**

```

  

2) typing the same command followed by another supported command as an additional aggument will result in performiing the second command without explicitly entering the container shell. For example:

```

[user@cn3316 ~]$ **ag python /ANNOgesic/tests/test\_gen\_svg.py** 
.
----------------------------------------------------------------------
Ran 1 test in 0.004s

OK
[user@cn3316 ~]$ **ag python3 /ANNOgesic/tests/test\_converter.py** 
..........
----------------------------------------------------------------------
Ran 10 tests in 0.018s

OK

```

In particular, the following command will display ANNOgesic help message:

```

[user@cn3316 ~]$ **ag annogesic --help** 

       ___    _   ___   ______                  _     
      /   |  / | / / | / / __ \____ ____  _____(_)____ \
  __ / /| | /  |/ /  |/ / / / / __ `/ _ \/ ___/ / ___/__\
 |  / ___ |/ /|  / /|  / /_/ / /_/ /  __(__  ) / /__    /
 | /_/  |_/_/ |_/_/ |_/\____/\__, /\___/____/_/\___/   /
 |                          /____/ 
 |__________________
 |_____________________
 |________________________________________________
 |                                                \
 |________________________________________________/

usage: annogesic [-h] [--version]
                 {create,get_input_files,update_genome_fasta,annotation_transfer,tss_ps,optimize_tss_ps,terminator,transcript,utr,srna,sorf,promoter,operon,circrna,go_term,srna_target,snp,ppi_network,localization,riboswitch_thermometer,crispr,merge_features,screenshot,colorize_screenshot_tracks}
                 ...

positional arguments:
  {create,get_input_files,update_genome_fasta,annotation_transfer,tss_ps,optimize_tss_ps,terminator,transcript,utr,srna,sorf,promoter,operon,circrna,go_term,srna_target,snp,ppi_network,localization,riboswitch_thermometer,crispr,merge_features,screenshot,colorize_screenshot_tracks}
                        commands
    create              Create a project
    get_input_files     Get required files. (i.e. annotation files, fasta
                        files)
    update_genome_fasta
                        Get fasta files of reference genomes if the reference
                        sequences do not exist.
    annotation_transfer
                        Transfer the annotations from a closely related
                        species genome to a target genome.
    tss_ps              Detect TSSs or processing sites.
    optimize_tss_ps     Optimize TSSs or processing sites based on manual
                        detected ones.
    terminator          Detect rho-independent terminators.
    transcript          Detect transcripts based on coverage file.
    utr                 Detect 5'UTRs and 3'UTRs.
    srna                Detect intergenic, antisense and UTR-derived sRNAs.
    sorf                Detect expressed sORFs.
    promoter            Discover promoter motifs.
    operon              Detect operons and sub-operons.
    circrna             Detect circular RNAs.
    go_term             Extract GO terms from Uniprot.
    srna_target         Detect sRNA-mRNA interactions.
    snp                 Detect SNP/mutation and generate fasta file if
                        mutations were found.
    ppi_network         Detect protein-protein interactions suported by
                        literature.
    localization        Predict subcellular localization of proteins.
    riboswitch_thermometer
                        Predict riboswitches and RNA thermometers.
    crispr              Predict CRISPR related RNAs.
    merge_features      Merge all features to one gff file.
    screenshot          Generate screenshots for selected features using IGV.
    colorize_screenshot_tracks
                        Add color information to screenshots (e.g. useful for
                        dRNA-Seq based TSS and PS detection. It only works
                        after running "screenshot" (after running batch
                        script).

optional arguments:
  -h, --help            show this help message and exit
  --version, -v         show version


```

Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. ANNOgesic.sh). For example:



```

#!/bin/bash
module load ANNOgesic
ag python  /ANNOgesic/tests/mock_gff3.py
ag python  /ANNOgesic/tests/mock_helper.py
ag python  /ANNOgesic/tests/test_blast_class.py
ag python  /ANNOgesic/tests/test_change_db_format.py
ag python  /ANNOgesic/tests/test_check_orphan.py
ag python3 /ANNOgesic/tests/test_circRNA.py
ag python  /ANNOgesic/tests/test_circrna.py
ag python  /ANNOgesic/tests/test_color_png.py
ag python  /ANNOgesic/tests/test_combine_frag_tex.py
ag python3 /ANNOgesic/tests/test_combine_gff.py
ag python  /ANNOgesic/tests/test_compare_sRNA_sORF.py
ag python3 /ANNOgesic/tests/test_converter.py

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] ANNOgesic.sh
```





