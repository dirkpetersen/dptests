

document.querySelector('title').textContent = 'randfold on Biowulf';
randfold on Biowulf


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




Randfold computes the probability that, for a given RNA sequence, the
Minimum Free Energy (MFE) of the secondary structure is different from a
distribution of MFE computed with random sequences obtained by permuting
the input sequence.



### References:


* Eric Bonnet, J. Wuyts, P. Rouzé and Y. Van de Peer. *Evidence that 
 microRNA precursors, unlike other non-coding RNAs, have lower folding
 free energies than random sequences.* Bioinformatics, 2004, 
 22:2911-2917.
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/15217813)
  | 
 PMC | 
 [Journal](http://bioinformatics.oxfordjournals.org/content/20/17/2911.long)



Documentation
* [Source 
 and supplemental data](http://bioinformatics.psb.ugent.be/software/details/Randfold)


Important Notes
* Module Name: randfold (see [the modules page](/apps/modules.html) for more information)



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int).



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load randfold**

```


Example: running randfold on a single miRNA sequence




```

[user@cn3144]$ **randfold**
FATAL: Usage: randfold <method> <file name> <number of randomizations>

methods:
-s simple mononucleotide shuffling
-d dinucleotide shuffling
-m markov chain 1 shuffling

Example: randfold -d let7.tfa 999

[user@cn3144]$ **cat > cel-let7.fa <<EOF**
>cel-let-7 Caenorhabditis elegans let-7 precursor RNA
UACACUGUGGAUCCGGUGAGGUAGUAGGUUGUAUAGUUUGGAAUAUUACCACCGGUGAAC
UAUGCAAUUUUCUACCUUACCGGAGACAGAACUCUUCGA
**EOF**
[user@cn3144]$ **randfold -d cel-let7.fa 999**
cel-let-7       -42.90  0.001000

```


Running randfold on a set of mouse miRNA sequences in series (randfold
will process one sequence at a time).




```

[user@cn3144]$ **wget "ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz"**
[user@cn3144]$ **gunzip hairpin.fa.gz**
[user@cn3144]$ **awk '/^>/ {p=0} /^>mmu/{p=1} p==1' hairpin.fa > mouse\_hairpin.fa**
[user@cn3144]$ **randfold -d mouse\_hairpin.fa 99 > mouse\_hairpin.randfold**

```

End the interactive session



```

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. randfold.sh), which uses the input file 'randfold.in'. For example:



```

#! /bin/bash
#SBATCH --job-name=randfold
set -e

module load randfold
inf=/data/$USER/test_data/randfold/mouse_hairpins.fa
outf=/data/$USER/test_data/randfold/mouse_hairpins.fa.randfold
randfold -d $inf 99 > $outf

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch randfold.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. randfold.swarm). For example:



```

randfold -d mouse_hairpins/mmu-let-7a-1.fasta 999 > mouse_hairpins/mmu-let-7a-1.fasta.randfold
randfold -d mouse_hairpins/mmu-let-7a-2.fasta 999 > mouse_hairpins/mmu-let-7a-2.fasta.randfold
randfold -d mouse_hairpins/mmu-let-7b.fasta 999 > mouse_hairpins/mmu-let-7b.fasta.randfold

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f randfold.swarm --module randfold
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module randfold  Loads the randfold module for each subjob in the swarm 
 | |
 | |
 | |








