

document.querySelector('title').textContent = 'OpenCRAVAT: a platform for the annotation of human genetic variation';
**OpenCRAVAT: a platform for the annotation of human genetic variation**


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



OpenCRAVAT is a new open source, scalable decision support system for
variant and gene prioritization. It includses a modular resource catalog 
to maximize community and developer involvement, and as a result the catalog
is being actively developed and growing every month. Resources made available via
the store are well-suited for analysis of cancer, as well as Mendelian and complex
diseases.



OB
### References:


* Kymberleigh A. Pagel, Rick Kim, Kyle Moad, Ben Busby, Lily Zheng, Collin Tokheim,
Michael Ryan, and Rachel Karchin   

*Integrated Informatics Analysis of Cancer-Related Variants*
 [JCO Clinical Cancer Informatics](https://ascopubs.org/doi/full/10.1200/CCI.19.00132) , 2020, v.4, p.310-317.


Documentation
* [OpenCRAVAT quick start guide](https://github.com/KarchinLab/open-cravat/wiki/quickstart-command-line)
* [OpenCRAVAT Github wiki page](https://github.com/KarchinLab/open-cravat/wiki)
* [OpenCRAVAT Home page](https://bio.tools/OpenCRAVAT)
* [Instructions for the GUI usage](https://github.com/KarchinLab/open-cravat/wiki/5.-GUI-usage)
* [OpenCRAVAT NCI Webinar 9/2021](https://www.youtube.com/watch?v=5Y4U1ocMs4U&ab_channel=OpenCRAVAT)


Important Notes
* Module Name: OpenCRAVAT (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **OC\_HOME**  OpenCRAVAT installation directory
	+ **OC\_BIN**       OpenCRAVAT executable directory
	+ **OC\_SRC**     OpenCRAVAT source code directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=16g --gres=gpu:p100:1,lscratch:10 -c4**
[user@cn2389 ~]$ **module load OpenCRAVAT** 
[+] Loading annovar 2020-06-08 on cn2389
[+] Loading OpenCRAVAT 2.3.1  ...

```

In order to annotate and interpret variants, OpenCRAVAT (OC) makes use of a database comprising chanks of data called "modules" (not to be confused with the Biowulf modules). The currenly installed modules OC modules are available in the folder: $OC\_MODULES.   
  

To get started with OpenCRAVAT, create an example input file in your current working directory by using the commands:

```

[user@cn2389 ~]$ **cd /data/$USER**
[user@cn2389 ~]$ **mkdir OpenCRAVAT && cd OpenCRAVAT**
[user@cn2389 ~]$ **oc new example-input .** 

```

The latter command will create a file "example\_input" in your current directory:

```

[user@cn2389 ~]$ **cat example\_input | wc -l** 
373
[user@cn2389 ~]$  **head -n 20 example\_input**
chr10   121593817       -       A       T       s0
chr10   2987654 +       T       A       s1
chr10   43077259        +       A       T       s2
chr10   8055656 +       A       T       s3
chr10   87864470        +       A       T       s4
chr10   87864486        +       A       -       s0
chr10   87864486        +       AA      -       s1
chr10   87894027        +       -       CG      s2
chr10   87894027        +       -       CT      s3
chr1    100719861       +       A       T       s4
chr1    10100   +       C       T       s0
chr1    110340653       +       CGGCTTT -       s1
chr11   108227625       +       A       T       s2
chr11   113789394       +       G       A       s3
chr1    111762684       +       G       A       s4
chr11   119206418       +       A       T       s0
chr1    114713881       +       TGGTC   -       s1
chr1    114713881       +       TGGTCTC -       s2
chr1    114716160       -       A       T       s3
chr11   1584916 +       -       GCC     s4

```

Run OpenCRAVAT on the test input file:

```

[user@cn2389 ~]$ **oc run ./example\_input -l hg38**
Input file(s): /vf/users/user/OpenCRAVAT/example_input
Genome assembly: hg38
Running converter...
        Converter (converter)           finished in 2.055s
Running gene mapper...                  finished in 21.534s
Running annotators...
        annotator(s) finished in 1.412s
Running aggregator...
        Variants                        finished in 0.174s
        Genes                           finished in 0.101s
        Samples                         finished in 0.192s
        Tags                            finished in 0.295s
Indexing
        variant base__chrom     finished in 0.009s
        variant base__coding    finished in 0.010s
        variant base__so        finished in 0.009s
Running postaggregators...
        Tag Sampler (tagsampler)        finished in 0.122s
Finished normally. Runtime: 26.077s

```

Once the job is finished, the following files wil be created:   
  


```

example_input.log  

example_input.sqlite  

example_input.err  


```

In particular, file example\_input.sqlite is the sqlite database with the results.
This sqlite database can be opened in the OpenCRAVAT web viewer as follows:

```

[user@cn2389 ~]$ **oc gui example\_input.sqlite**

   ____                   __________  ___ _    _____  ______
  / __ \____  ___  ____  / ____/ __ \/   | |  / /   |/_  __/
 / / / / __ \/ _ \/ __ \/ /   / /_/ / /| | | / / /| | / /
/ /_/ / /_/ /  __/ / / / /___/ _, _/ ___ | |/ / ___ |/ /
\____/ .___/\___/_/ /_/\____/_/ |_/_/  |_|___/_/  |_/_/
    /_/

OpenCRAVAT is served at 0.0.0.0:8080
(To quit: Press Ctrl-C or Ctrl-Break if run on a Terminal or Windows, or click "Cancel" and then "Quit" if run through OpenCRAVAT app on Mac OS)
START /usr/bin/exo-open --launch WebBrowser "http://0.0.0.0:8080/result/index.html?dbpath=/vf/users/user/OpenCRAVAT/example_input.sqlite"
...

```

On your local system, open a new window and type:

```

**ssh -t -L 8080:localhost:8080 biowulf.nih.gov "ssh -L 8080:localhost:8080 cn2389"**

```

Navigate a browser on your local system to the URL: localhost:8080.   
  

![](opencravat/opencravat_gui.png)
  
   

Exit an interactive session:

```

[user@cn2389 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```








