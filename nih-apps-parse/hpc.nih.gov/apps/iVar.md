

document.querySelector('title').textContent = 'iVar: iVar: a tool for viral amplicon-based sequencing.';
iVar: iVar: a tool for viral amplicon-based sequencing.


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



iVar is a computational package that contains functions broadly useful for viral amplicon-based sequencing. Additional tools for metagenomic sequencing are actively being incorporated into iVar. While each of these functions can be accomplished using existing tools, iVar contains an intersection of functionality from multiple tools that are required to call iSNVs and consensus sequences from viral sequencing data across multiple replicates.



### References:


* Nathan D. Grubaugh, Karthik Gangavarapu, Joshua Quick, Nathaniel L. Matteson, Jaqueline Goes De Jesus, Bradley J. Main, Amanda L. Tan, Lauren M. Paul, Doug E. Brackney, Saran Grewal, Nikos Gurfield, Koen K. A. Van Rompay, Sharon Isern, Scott F. Michael, Lark L. Coffey, Nicholas J. Loman and Kristian G. Andersen   

 *An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar.*  

[Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1618-7)  (2019), **5**(8). https://doi.org/10.1186/s13059-018-1618-7


Documentation
* [iVar Github page](https://github.com/andersen-lab/ivar)
* [iVar manual](https://andersen-lab.github.io/ivar/html/manualpage.html)


Important Notes
* Module Name: iVar (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **IVAR\_HOME**  IVAR installation directory
	+ **IVAR\_BIN**       IVAR executable directory
	+ **IVAR\_SRC**     a folder containing the source code
	+ **IVAR\_DATA**    a folder containing sample data
	+ **IVAR\_HMM**     a folder containing sample Hidden Markov Models



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive** 
[user@cn4471 ~]$ **module load ivar** 
[+] Loading samtools 1.9  ...
[+] Loading ivar 1.3.1  ...

```

Copy sample data to your current folder:

```

[user@cn4471 ~]$ **cp -r $IVAR\_DATA/\* .**
[user@cn4471 ~]$ **cp test\_amplicon.sorted.bam test.bam** 

```

Preprocess the data with samtools:
[user@cn4471 ~]$ **samtools sort -o test.sorted.bam test.bam** 
[user@cn4471 ~]$ **samtools index test.sorted.bam**

iVar can be used as follows:

```

[user@cn4471 ~]$ **ivar -h**
Usage:  ivar [command <trim|callvariants|filtervariants|consensus|getmasked|removereads|version|help>]

        Command       Description
           trim       Trim reads in aligned BAM file
       variants       Call variants from aligned BAM file
 filtervariants       Filter variants across replicates
      consensus       Call consensus from aligned BAM file
      getmasked       Detect primer mismatches and get primer indices for the amplicon to be masked
    removereads       Remove reads from trimmed BAM file
        version       Show version information

To view detailed usage for each command, type **ivar <command>**, for example:
[user@cn4471 ~]$ **ivar trim -h** 
Usage: ivar trim -i <input.bam> -b <primers.bed> -p <prefix> [-m <min-length>] [-q <min-quality>] [-s <sliding-window-width>]
Input Options    Description
           -i    (Required) Sorted bam file, with aligned reads, to trim primers and quality
           -b    (Required) BED file with primer sequences and positions
           -m    Minimum length of read to retain after trimming (Default: 30)
           -q    Minimum quality threshold for sliding window to pass (Default: 20)
           -s    Width of sliding window (Default: 4)
           -e    Include reads with no primers. By default, reads with no primers are excluded
Output Options   Description
           -p    (Required) Prefix for the output BAM file

```

Perform the trimming on sample data:

```

[user@cn4471 ~]$ **cp test\_isize.bed test\_primers.bed** 
[user@cn4471 ~]$ **ivar trim -b test\_primers.bed -p test.trimmed -i test.bam -q 15 -m 50 -s 4** 
Found 218 primers in BED file
Amplicons detected:

Number of references in file: 1
NC_045512.2
Using Region: NC_045512.2

Found 8 mapped reads
Found 0 unmapped reads
Sorted By Coordinate
-------
Processed 10% reads ...
Processed 20% reads ...
Processed 30% reads ...
Processed 40% reads ...

-------
Results:
Primer Name     Read Count
nCoV-2019_1_LEFT        0
nCoV-2019_1_RIGHT       0
nCoV-2019_2_LEFT        0
nCoV-2019_2_RIGHT       0
nCoV-2019_3_LEFT        0
nCoV-2019_3_RIGHT       0
nCoV-2019_4_LEFT        0
nCoV-2019_4_RIGHT       0
...
nCoV-2019_96_LEFT       0
nCoV-2019_96_RIGHT      0
nCoV-2019_97_LEFT       0
nCoV-2019_97_RIGHT      0
nCoV-2019_98_LEFT       0
nCoV-2019_98_RIGHT      0

Trimmed primers from 100% (8) of reads.
0% (0) of reads were quality trimmed below the minimum length of 50 bp and were not written to file.
0% (0) of reads that started outside of primer regions were not written to file
50% (4) of reads had their insert size smaller than their read length

```

End the interactive session:

```

[user@cn4471 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





