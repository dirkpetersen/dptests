

document.querySelector('title').textContent = 'scChromHMM: a suite of tools for rapid processing of single-cell histone modification data. ';
**scChromHMM: a suite of tools for rapid processing of single-cell histone modification data.** 


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



scChromHMM provides a suite of tools for rapid processing of single-cell histone modification data to perform chromatin states analysis of the genome within each single-cell. It is an extention of bulk ChromHMM framework, which consumes the HMM model learned from ChromHMM and perform chromatin state analysis by running forward-backward algorithm for each single-cell.



### References:


* Bingjie Zhang, Avi Srivastava, Eleni Mimitou, Tim Stuart, Ivan Raimondi, Yuhan Hao, Peter Smibert & Rahul Satija   

*Characterizing cellular heterogeneity in chromatin state with scCUT&Tag-pro*   

[Nature Biotechnology, 24 March 2022](https://www.nature.com/articles/s41587-022-01250-0)
Documentation+ [scChromHMM Github page](https://github.com/satijalab/scChromHMM)

Important Notes+ Module Name: scChromHMM (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
+ Unusual environment variables set
	- **SCCHROMHMM\_HOME**  installation directory


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
  


```

[user@biowulf]$ **sinteractive** 
[user@cn0861 ~]$ **module scChromHMM** 
[+] Loading singularity  3.8.5-1  on cn4177
[+] Loading scChromHMM 20220610  ...

```


```

[user@cn0861 ~]$ **schrom -h** 
chrom 0.1.0
Avi Srivastava, Bingjie Zhang, Rahul Satija
Generate summary stats for multimodal data.

USAGE:
    schrom [SUBCOMMAND]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

SUBCOMMANDS:
    help         Prints this message or the help of the given subcommand(s)
    hmm          A subcommand to run hmm.
    transform    A subcommand to transform long form matrices to short.
[user@cn0861 ~]$ **schrom hmm -h** 
schrom-hmm
A subcommand to run hmm.

USAGE:
    schrom hmm [FLAGS] --anchors ... --common\_cells  --fragments ... --model  --output  --threads 

FLAGS:
 -h, --help Prints help information
 --onlyone quantify only chromosome one
 -V, --version Prints version information

OPTIONS:
 -a, --anchors ... path to the anchors files. [Same order as fragments]
 -c, --common\_cells  path to the file with cellular barcodes of common assay
 -f, --fragments ... path to the fragment files.
 -m, --model  path to the chromeHMM model.txt file
 -o, --output  path to the output directory
 -t, --threads  number of threads to use
[user@cn0861 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226

```





