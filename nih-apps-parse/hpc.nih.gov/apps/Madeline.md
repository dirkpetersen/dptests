

document.querySelector('title').textContent = 'Madeline: pedigree drawing with an emphasis on readability and aesthetics.';
**Madeline: pedigree drawing with an emphasis on readability and aesthetics.**


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



The Madeline 2.0 Pedigree Drawing Engine (PDE) is a
pedigree drawing program for use in linkage and family-based
association studies. The program is designed to handle large and
complex pedigrees with an emphasis on readability and aesthetics.



### References:


* Edward H. Trager, Ritu Khanna, Adrian Marrs, Lawrence Siden,
Kari E.H. Branham, Anand Swaroop and Julia E. Richards,   

*Madeline 2.0 PDE: a new program for local and web-based pedigree drawing.* [Bioinformatics](https://academic.oup.com/bioinformatics/article/23/14/1854/189087?login=true), 2007, 
**23**(14), pp. 1854–1856.


Documentation
* [Madeline 2.0\_PDE Github page](https://github.com/piratical/Madeline_2.0_PDE)
* [Madeline 2.0\_PDE Home page](https://madeline.med.umich.edu/madeline/index.php)


Important Notes
* Module Name: Madeline (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **MADELINE\_HOME**  MADELINE installation directory
	+ **MADELINE\_BIN**       MADELINE executable directory
	+ **MADELINE\_SRC**     MADELINEsource code directory
	+ **MADELINE\_DATA**     MADELINE sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive** 
[user@cn0849 ~]$ **module load Madeline** 
[+] Loading madeline  2.0  on cn0849
[+] Loading singularity  3.8.5-1  on cn0849
[+] Loading ImageMagick  7.0.7  on cn0849

```

Copy test data to your current filder:

```

[user@cn0849 ~]$ **cp -r $MADELINE\_DATA/\* .** 

```

Run the madeline executable on the test data:

```

[user@cn0849 ~]$ **madeline -L "IndividualID DOB" ./input/si\_001.data**
┌─────────────────────────────┐
│ Welcome to Madeline 2.0 PDE │
└─────────────────────────────┘
--------------------------------------------
 LABELS                          TOTAL: 2
--------------------------------------------
 1. IndividualID
 2. DOB
--------------------------------------------
Parser::readFile(): Opening a(n) UTF-8 file ...
Reading file data in Madeline flat file format ...
----------------------------------------------------------------
 ICON COLUMNS                                        TOTAL: 1
----------------------------------------------------------------
   1. Affected                         has 2 non-missing levels.
----------------------------------------------------------------
Table 1 is a pedigree table.
  Start of    addPedigreesFromDataTable 
Warning: Adding virtual father “S00100” who is not present in the input data file.
Warning: Adding virtual mother “S00101” who is not present in the input data file.
Siblings are ordered by DOB.
  End of      addPedigreesFromDataTable 
  Start of    draw                       
Pedigree output file is “si_001_pedigree.svg”
  End of      draw                       

[user@cn0849 ~]$ **madeline --color ./input/si\_002.data**
┌─────────────────────────────┐
│ Welcome to Madeline 2.0 PDE │
└─────────────────────────────┘
--------------------------------------------
 LABELS                          TOTAL: 0
--------------------------------------------
--------------------------------------------
Parser::readFile(): Opening a(n) UTF-8 file ...
Reading file data in Madeline flat file format ...
----------------------------------------------------------------
 ICON COLUMNS                                        TOTAL: 1
----------------------------------------------------------------
   1. Affected                         has 3 non-missing levels.
----------------------------------------------------------------
Table 1 is a pedigree table.
  Start of    addPedigreesFromDataTable  
Siblings are ordered by DOB.
  End of      addPedigreesFromDataTable  
  Start of    draw                       
Pedigree output file is “si_002_pedigree.svg”
  End of      draw                       

```

Visualize the results: 

```

[user@cn0849 ~]$ **display si\_001\_pedigree.svg**

```

  

![](madeline/si_001_pedigree.svg)
  
  
  


```

[user@cn0849 ~]$ **display si\_002\_pedigree.svg**

```

![](madeline/si_002_pedigree.svg) 
[user@cn0849 ~]$ **exit**
exit
salloc.exe: Relinquishing job allocation 28896415
salloc.exe: Job allocation 28896415 has been revoked.
[user@biowulf ~]$








