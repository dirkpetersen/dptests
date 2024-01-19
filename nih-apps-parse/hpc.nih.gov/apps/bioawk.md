

document.querySelector('title').textContent = 'bioawk on Biowulf &amp; Helix';
bioawk on Biowulf & Helix



|  |
| --- |
| 
Quick Links
[Interactive job on Biowulf](#int)
[Batch job on Biowulf](#serial)
[Swarm of jobs](#swarm)
 |



Description
Bioawk extends awk with support for several common biological data formats,
including optionally gzip'ed BED, GFF, SAM, VCF, FASTA/Q and TAB-delimited
formats with column names. It also adds a few built-in functions and an command
line option to use TAB as the input/output delimiter. When the new
functionality is not used, bioawk is intended to behave exactly the same as the
original BWK awk. 


There may be multiple versions of bioawk available. An easy way of selecting the
version is to use [modules](/apps/modules.html). To see the modules
available, type



```

module avail bioawk 

```

To select a module use



```

module load bioawk/[version]

```

where `[version]` is the version of choice.


### Environment variables set


* `$PATH`
* `$BIOAWK_TEST_DATA`


### Documentation


* [GitHub page](https://github.com/lh3/bioawk)




Interactive job on Biowulf
Allocate an interactive session with [sinteractive](/docs/userguide.html#int)
and use as described below



```

biowulf$ **sinteractive** 
node$ **module load bioawk samtools**
node$ # what formats are supported and what automatic variables are created for each of them?
node$ **bioawk -c help**
bed:
        1:chrom 2:start 3:end 4:name 5:score 6:strand 7:thickstart 8:thickend 9:rgb 10:blockcount 11:blocksizes 12:blockstarts 
sam:
        1:qname 2:flag 3:rname 4:pos 5:mapq 6:cigar 7:rnext 8:pnext 9:tlen 10:seq 11:qual 
vcf:
        1:chrom 2:pos 3:id 4:ref 5:alt 6:qual 7:filter 8:info 
gff:
        1:seqname 2:source 3:feature 4:start 5:end 6:score 7:filter 8:strand 9:group 10:attribute 
fastx:
        1:name 2:seq 3:qual 4:comment 
node$ # create fasta from bam file. Reverse complement the read sequence if the
node$ # read aligned to the minus strand
node$ **samtools view $BIOAWK\_TEST\_DATA/aln.bam \
 | bioawk -c sam '{s=$seq; if(and($flag, 16)) {s=revcomp($seq)} print ">"$qname"\n"s}' \
 | head**
>D00777:83:C8R9PACXX:3:2111:6151:30663
CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAAACCCTAACACTGACGCGAACGGAAAGAGTAA
>D00777:83:C8R9PACXX:3:2111:6151:30663
CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAAACCCTAACACTGACGCGAACGGAAAGAGTAA
>D00777:83:C8R9PACXX:3:2111:6151:30663
CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAAACCCTAACACTGACGCGAACGGAAAGAGTAA
>D00777:83:C8R9PACXX:3:2111:6151:30663
CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAAACCCTAACACTGACGCGAACGGAAAGAGTAA
>D00777:83:C8R9PACXX:3:2111:6151:30663
CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAAACCCTAACACTGACGCGAACGGAAAGAGTAA
node$ **samtools view $BIOAWK\_TEST\_DATA/aln.bam | bioawk -c sam '$seq ~ /^ATG/ {printf("%s:%i\t%s\n", $rname, $pos, $seq)}'**
chr1:14686      ATGGAGCCCCCTACGATTCCCAGTCGTCCTCGTCCTCCTCTGCCTGTGGCTGCTGCGGTGGCGGCAGAGGAGGGATGGAGTCTGACACGCGGGCAAAGGCT
chr1:17535      ATGCCCTGGGTCCCCACTAAGCCAGGCCGGGCCTCCCGCCCACACCCCTCGGCCCTGCCCTCTGGCCATACAGGTTCTCGGTGGTGTTGAAGAGCAGCAAG
chr1:17535      ATGCCCTGGGTCCCCACTAAGCCAGGCCGGGCCTCCCGCCCACACCCCTCGGCCCTGCCCTCTGGCCATACAGGTTCTCGGTGGTGTTGAAGAGCAGCAAG
node$ **exit**
biowulf$

```



Batch job on Biowulf
Create a batch script similar to the following example:



```

#! /bin/bash
# this file is bioawk.batch

module load bioawk || exit 1
bioawk -c gff '$feature == "exon" and (end - start) < 100' refseq.gff

```

Submit to the queue with [sbatch](/docs/userguide.html):



```

biowulf$ **sbatch --time=10 bioawk.batch**

```



Swarm of jobs on Biowulf
Create a swarm command file similar to the following example:



```

# this file is bioawk.swarm
bioawk -c fastx '{print ">" $name; print revcomp($seq)}' seq1.fa.gz | gzip -c > seq1.rc.fa.gz
bioawk -c fastx '{print ">" $name; print revcomp($seq)}' seq2.fa.gz | gzip -c > seq2.rc.fa.gz
bioawk -c fastx '{print ">" $name; print revcomp($seq)}' seq3.fa.gz | gzip -c > seq3.rc.fa.gz

```

And submit to the queue with [swarm](/apps/swarm.html)



```

biowulf$ **swarm -f bioawk.swarm --module bioawk --time 10**

```



