

document.querySelector('title').textContent = 'netOglyc on Biowulf';
netOglyc on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int)
[Batch job](#sbatch)
[Swarm of jobs](#swarm)
[Multi-sequence fasta file](#multi)
 |



NetOglyc produces neural network predictions of mucin type GalNAc O-glycosylation sites in mammalian proteins.



### References:


* K. Julenius, A. Mølgaard, R. Gupta and S. Brunak.
 [Prediction, conservation analysis and structural characterization of mammalian mucin-type O-glycosylation sites](https://www.ncbi.nlm.nih.gov/pubmed/?term=15385431)
*Glycobiology,* 15:153-164, 2005.


Documentation
* [netOglyc 4.0 Main Site](https://services.healthtech.dtu.dk/services/NetOGlyc-4.0/)


Important Notes
* Module Name: netOglyc (see [the modules page](/apps/modules.html) for more information)
* singlethreaded
 * Environment variables set 
	+ NETOGLYC\_HOME
	+ NETOGLYC\_EXAMPLES* Example files in $NETOGLYC\_EXAMPLES
* **This application requires /lscratch be allocated (see below)**
* The input fasta file header must conform to the following rules:
	+ No spaces are allowed between the inital > and the description:  
	
	**GOOD:** 
	```
	>NM_001008540.2 Homo sapiens C-X-C motif
	```
	
	**BAD:** 
	```
	>   NM_001008540.2 Homo sapiens C-X-C motif
	```
	+ The vertical bar character | is not allowed in the description:  
	
	**GOOD:** 
	```
	>NM_000758.4 Homo sapiens colony stimulating factor 2 (CSF2), mRNA
	```
	
	**BAD:** 
	```
	>NM_000758.4 Homo sapiens colony stimulating factor 2 | other stuff
	```
* While netOglyc can process a multiple sequence fasta file, there is a limit to the input size of the file. It is better to split the fasta file into single sequence files prior to running. See below for an example perl script for doing this on the fly.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --gres=lscratch:10**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load netOglyc**
[user@cn3144 ~]$ **netOglyc $NETOGLYC\_EXAMPLES/GLP\_MACFU.fsa**
##gff-version 2
##source-version NetOGlyc 4.0.0.11
##date 23-9-6
##Type Protein
#seqname        source  feature start   end     score   strand  frame   comment
GLP_MACFU       netOGlyc-4.0.0.11       CARBOHYD        1       1       0.680709        .       .       #POSITIVE
GLP_MACFU       netOGlyc-4.0.0.11       CARBOHYD        2       2       0.790723        .       .       #POSITIVE
GLP_MACFU       netOGlyc-4.0.0.11       CARBOHYD        3       3       0.848504        .       .       #POSITIVE
GLP_MACFU       netOGlyc-4.0.0.11       CARBOHYD        4       4       0.707939        .       .       #POSITIVE...
...

```

If something goes wrong, temporary files (including logs) are written to and compressed in /lscratch/$SLURM\_JOB\_ID:

```

[user@cn3144 ~]$ **ls /lscratch/$SLURM\_JOB\_ID**
netOGlyc-1474660.tar.gz  netOGlyc-1474970.tar.gz

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. netOglyc.sh). For example:



```

#!/bin/bash
set -e
module load netOglyc
netOglyc my_fasta_file.fasta > my_fasta_file.out

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --gres=lscratch:10 netOglyc.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. netOglyc.swarm). For example:



```

netOglyc < 1.fasta > 1.out
netOglyc < 2.fasta > 2.out
netOglyc < 3.fasta > 3.out
netOglyc < 4.fasta > 4.out

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f netOglyc.swarm --module netOglyc --gres lscratch:10
```

where


|  |  |
| --- | --- |
| --module netOglyc Loads the netOglyc module for each subjob in the swarm  | |
| --gres lscratch:10 Allocates 10 GB of /lscratch | |


Multi-sequence fasta break-up script
Very large, multi-sequence fasta files are not handled well with netOglyc. It is safer to keep the fasta input small.


Here is a perl script that will take a multi-sequence fasta file, break it up into single fasta files, and run netOglyc on each. The output is appended to a single output file.



```

#!/usr/bin/perl

use File::Temp qw/ tempfile /;
my (undef, $infile) = tempfile();

my $ii=0;
my $i=0;
my $order="";
my $tim=0;
my $com="netOglyc $infile.fsa >& $infile.out";
my $nam="";
my $seq="";
my $INPUT = $ARGV[0];
my $OUTPUT = $ARGV[1];

open (CFGFILE, $ARGV[0]);
while (){
 chomp;
 my $line = $\_;
 if ($line =~ />/) {
 $i++;
 $seq[$i]=$line."\n";
 $nam[$i]=$line;
 } else {
 $seq[$i] = $seq[$i].$line."\n";
 }
}

unlink $OUTPUT;

while ($ii<$i){
 $ii++;
 open OUTFILE, "> $infile.fsa";
 $order=$seq[$ii];
 print "Sent ".$order." ".$ii." of ".$i." to NetOglyc"."\n"."\n";
 print OUTFILE $order;
 close OUTFILE;
 open OUTFILE, ">> $ARGV[1]";
 $order=$nam[$ii];
 print OUTFILE $order."\n";
 close OUTFILE;
 unlink "$infile.out";
 system ($com);

 open INFILE, "$infile.out";
 open OUTFILE, ">> $ARGV[1]";
 while (){
 print OUTFILE $\_;
 }
 print OUTFILE "\n";
 close OUTFILE;
}
unlink "$infile.out";
unlink "$infile.fsa";

```

The script can be run like so:



```
perl multi.pl multi-fasta.fsa multi-fasta.out
```









