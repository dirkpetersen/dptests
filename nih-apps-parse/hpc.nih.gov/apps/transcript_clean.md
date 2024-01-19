

document.querySelector('title').textContent = "TEMPLATE";
Transcript Clean on Biowulf


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



Description



### References:


* Blow J.
 [**A really amazing research paper.**](https://www.ncbi.nlm.nih.gov/pubmed/00000000)
*J Mol Biol. 2012 Jan 13;415(2):406-18.*
* Blow J., Doe J.
 [**A retread of another amazing research paper.**](https://www.ncbi.nlm.nih.gov/pubmed/00000000)
*J Struct Biol. 2012 Dec;180(3):519-30.*


Documentation
* [TEMPLATE Main Site](https://hpcwebdev.cit.nih.gov/)


Important Notes
This application requires a [graphical connection using NX](/docs/connect.html#nx)


* Module Name: TEMPLATE (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded/singlethreaded/MPI...
 * This application produces HTML reports. You can use [hpcdrive](/docs/hpcdrive.html) to view these reports on your local workstation.
* Environment variables set 
	+ TEMPLATE\_HOME* Example files in ???* Reference data in /fdb/TEMPLATE/



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cncn4338 ~]$ **module load transcript\_clean**
[+] Loading transcript_clean  2.0.3  on cn4338 
[+] Loading singularity  3.10.5  on cn4338 
[user@cn4338 ~]$  **cd /usr/local/apps/transcript\_clean/2.0.3/TranscriptClean-2.0.3**
[user@cn4338 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226

```


Example
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Example using help command:



```

[user@cn4338 TranscriptClean-2.0.3]$ **python TranscriptClean.py --help** 
Usage: TranscriptClean.py [options]

Options:
  -h, --help            show this help message and exit
  -s FILE, --sam=FILE   Input SAM file containing transcripts to correct. Must
                        contain a header.
  -g FILE, --genome=FILE
                        Reference genome fasta file. Should be the same one
                        used during mapping to generate the provided SAM file.
  -t N_THREADS, --threads=N_THREADS
                        Number of threads to run program with.
  -j FILE, --spliceJns=FILE
                        Splice junction file obtained by mapping Illumina
                        reads to the genome using STAR, or alternately,
                        extracted from a GTF using the accessory script. More
                        formats may be supported in the future.
  -v FILE, --variants=FILE
                        VCF formatted file of variants to avoid correcting
                        away in the data (optional).
  --maxLenIndel=MAXLENINDEL
                        Maximum size indel to correct (Default: 5 bp)
  --maxSJOffset=MAXSJOFFSET
                        Maximum distance from annotated splice junction to
                        correct (Default: 5 bp)
  -o FILE, --outprefix=FILE
                        Output file prefix. '_clean' plus a file extension
                        will be added to the end.
  -m CORRECTMISMATCHES, --correctMismatches=CORRECTMISMATCHES
                        If set to false, TranscriptClean will skip mismatch
                        correction. Default: true
  -i CORRECTINDELS, --correctIndels=CORRECTINDELS
                        If set to false, TranscriptClean will skip indel
                        correction. Default: true
  --correctSJs=CORRECTSJS
                        If set to false, TranscriptClean will skip splice
                        junction correction. Default: true, but you must
                        provide a splice junction annotation file in order for
                        it to work.
  --dryRun              If this option is set, TranscriptClean will read in
                        the sam file and record all insertions, deletions, and
                        mismatches, but it will skip correction. This mode is
                        useful for checking the distribution of transcript
                        errors in the data before running correction.
  --primaryOnly         If this option is set, TranscriptClean will only
                        output primary mappings of transcripts (ie it will
                        filter                       out unmapped and
                        multimapped lines from the SAM input.
  --canonOnly           If this option is set, TranscriptClean will output
                        only canonical transcripts and transcripts containing
                        annotated noncanonical junctions to the clean SAM file
                        at the end of the run.
  --tmpDir=TMP_PATH     If you would like the tmp files to be written
                        somewhere different than the final output, provide the
                        path to that location here.
  --bufferSize=BUFFER_SIZE
                        Number of lines to output to file at once by each
                        thread during run. Default = 100
  --deleteTmp           If this option is set, the temporary directory
                        generated by TranscriptClean (TC_tmp) will be removed
                        at the end of the run.
    
```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. TEMPLATE.sh). For example:



```

#!/bin/bash
set -e
module load TEMPLATE
TEMPLATE < TEMPLATE.in > TEMPLATE.out

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] TEMPLATE.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. TEMPLATE.swarm). For example:



```

TEMPLATE < TEMPLATE.in > TEMPLATE.out
TEMPLATE < TEMPLATE.in > TEMPLATE.out
TEMPLATE < TEMPLATE.in > TEMPLATE.out
TEMPLATE < TEMPLATE.in > TEMPLATE.out

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f TEMPLATE.swarm [-g #] [-t #] --module TEMPLATE
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module TEMPLATE Loads the TEMPLATE module for each subjob in the swarm 
 | |
 | |
 | |










