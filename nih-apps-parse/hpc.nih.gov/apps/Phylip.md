

document.querySelector('title').textContent = 'Phylip on Biowulf';
Phylip on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



PHYLIP (PHYlogeny Inference Package) is a package of programs for inferring phylogenies (evolutionary trees). Methods that are available in the package include parsimony, distance matrix, and likelihood methods, including bootstrapping and consensus trees. Data types that can be handled include molecular sequences, gene frequencies, restriction sites and fragments, distance matrices, and discrete characters.



PHYLIP is a single-threaded program, and is intended to be used interactively on helix. The PHYLIP package consists of a large number of individual programs. They are installed in /usr/local/phylip/exe. [List of PHYLIP programs.](http://evolution.genetics.washington.edu/phylip/programs.html)



### References:


* [A primer to phylogenetic analysis using PHYLIP](http://koti.mbnet.fi/tuimala/oppaat/phylip2.pdf)


Documentation
* [Phylip Main Site](http://evolution.genetics.washington.edu/phylip.html)


Important Notes
* Module Name: phylip (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded (Interactive)
* Example files in /usr/local/apps/phylip/example



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **cp -r /usr/local/apps/phylip/example .**

[user@cn3144 ~]$ **cd example**

[user@cn3144 ~]$ **module load phylip**

[user@cn3144 ~]$ **dnapars**
dnapars: can't find input file "infile"
Please enter a new file name> seq

DNA parsimony algorithm, version 3.65

Setting for this run:
  U                 Search for best tree?  Yes
  S                        Search option?  More thorough search
  V              Number of trees to save?  10000
  J   Randomize input order of sequences?  No. Use input order
  O                        Outgroup root?  No, use as outgroup species  1
  T              Use Threshold parsimony?  No, use ordinary parsimony
  N           Use Transversion parsimony?  No, count all steps
  W                       Sites weighted?  No
  M           Analyze multiple data sets?  No
  I          Input sequences interleaved?  Yes
  0   Terminal type (IBM PC, ANSI, none)?  ANSI
  1    Print out the data at start of run  No
  2  Print indications of progress of run  Yes
  3                        Print out tree  Yes
  4          Print out steps in each site  No
  5  Print sequences at all nodes of tree  No
  6       Write out trees onto tree file?  Yes

  Y to accept these or type the letter for one to change
Y
Adding species:
   1. Archaeopt 
   2. Hesperorni
   3. Baluchithe
   4. B. virgini
   5. Brontosaur
   6. B.subtilis

Doing global rearrangements on all trees tied for best
  !-----------!
   ...........
   ...........

Collapsing best trees
Output written to file "outfile"
Tree also written onto file "outtree"
Done.

[user@cn3144 ~]$ **cat outfile**


DNA parsimony algorithm, version 3.65

One most parsimonious tree found:

                +---------------B.subtilis
  +-------------4  
  |             +----------Brontosaur
  |  
  |                      +-B. virgini
  |        +-------------3  
  1--------2             +------Baluchithe
  |        |  
  |        +-----------Hesperorni
  |  
  +-----------Archaeopt 

requires a total of     21.000

  between      and       length
  -------      ---       ------
     1           4       0.230769
     4      B.subtilis   0.269231
     4      Brontosaur   0.192308
     1           2       0.153846
     2           3       0.230769
     3      B. virgini   0.038462
     3      Baluchithe   0.115385
     2      Hesperorni   0.192308
     1      Archaeopt    0.192308

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. phylip.sh). For example:



```

#!/bin/bash
cd /data/$USER/somedir

module load phylip/3.696

dnaml << EOF
sequences.dat
Y
EOF

```


Basically, you need to provide Phylip with the same input that it would expect if you ran the program interactively. You can include the parameters directly in the batch script, as in the example above. Or you can put these parameters 
into a file (e.g. 'phylip\_input') and include:

```

dnaml < phylip_input > output

```

in your batch script. 



Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] phylip.sh
```







