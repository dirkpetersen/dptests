

document.querySelector('title').textContent = ' vcflib on Biowulf';

vcflib on Biowulf



|  |
| --- |
| 
Quick Links
[Interactive job](#int)
[Batch job on Biowulf](#serial)
[Parallel job on Biowulf](#parallel)
[Swarm of jobs](#swarm)
 |



Description
vcflib is a C++ library for parsing Variant Call Format (VCF) files and a
set of command line tools based on that library.


The following tools are currently available:



[vcf2dag](https://github.com/ekg/vcflib#vcf2dag), 
[vcf2fasta](https://github.com/ekg/vcflib#vcf2fasta), 
[vcf2tsv](https://github.com/ekg/vcflib#vcf2tsv), 
[vcfaddinfo](https://github.com/ekg/vcflib#vcfaddinfo), 
[vcfafpath](https://github.com/ekg/vcflib#vcfafpath), 
[vcfallelicprimitives](https://github.com/ekg/vcflib#vcfallelicprimitives), 
[vcfaltcount](https://github.com/ekg/vcflib#vcfaltcount), 
[vcfannotate](https://github.com/ekg/vcflib#vcfannotate), 
[vcfannotategenotypes](https://github.com/ekg/vcflib#vcfannotategenotypes), 
[vcfbreakmulti](https://github.com/ekg/vcflib#vcfbreakmulti), 
[vcfcat](https://github.com/ekg/vcflib#vcfcat), 
[vcfcheck](https://github.com/ekg/vcflib#vcfcheck), 
[vcfclassify](https://github.com/ekg/vcflib#vcfclassify), 
[vcfcleancomplex](https://github.com/ekg/vcflib#vcfcleancomplex), 
[vcfcombine](https://github.com/ekg/vcflib#vcfcombine), 
[vcfcommonsamples](https://github.com/ekg/vcflib#vcfcommonsamples), 
[vcfcountalleles](https://github.com/ekg/vcflib#vcfcountalleles), 
[vcfcreatemulti](https://github.com/ekg/vcflib#vcfcreatemulti), 
[vcfdistance](https://github.com/ekg/vcflib#vcfdistance), 
[vcfecho](https://github.com/ekg/vcflib#vcfecho), 
[vcfentropy](https://github.com/ekg/vcflib#vcfentropy), 
[vcfevenregions](https://github.com/ekg/vcflib#vcfevenregions), 
[vcffilter](https://github.com/ekg/vcflib#vcffilter), 
[vcffixup](https://github.com/ekg/vcflib#vcffixup), 
[vcfflatten](https://github.com/ekg/vcflib#vcfflatten), 
[vcfgeno2alleles](https://github.com/ekg/vcflib#vcfgeno2alleles), 
[vcfgeno2haplo](https://github.com/ekg/vcflib#vcfgeno2haplo), 
[vcfgenosamplenames](https://github.com/ekg/vcflib#vcfgenosamplenames), 
[vcfgenosummarize](https://github.com/ekg/vcflib#vcfgenosummarize), 
[vcfgenotypecompare](https://github.com/ekg/vcflib#vcfgenotypecompare), 
[vcfgenotypes](https://github.com/ekg/vcflib#vcfgenotypes), 
[vcfglbound](https://github.com/ekg/vcflib#vcfglbound), 
[vcfglxgt](https://github.com/ekg/vcflib#vcfglxgt), 
[vcfhetcount](https://github.com/ekg/vcflib#vcfhetcount), 
[vcfhethomratio](https://github.com/ekg/vcflib#vcfhethomratio), 
[vcfindex](https://github.com/ekg/vcflib#vcfindex), 
[vcfinfo2qual](https://github.com/ekg/vcflib#vcfinfo2qual), 
[vcfinfosummarize](https://github.com/ekg/vcflib#vcfinfosummarize), 
[vcfintersect](https://github.com/ekg/vcflib#vcfintersect), 
[vcfkeepgeno](https://github.com/ekg/vcflib#vcfkeepgeno), 
[vcfkeepinfo](https://github.com/ekg/vcflib#vcfkeepinfo), 
[vcfkeepsamples](https://github.com/ekg/vcflib#vcfkeepsamples), 
[vcfleftalign](https://github.com/ekg/vcflib#vcfleftalign), 
[vcflength](https://github.com/ekg/vcflib#vcflength), 
[vcfnumalt](https://github.com/ekg/vcflib#vcfnumalt), 
[vcfoverlay](https://github.com/ekg/vcflib#vcfoverlay), 
[vcfparsealts](https://github.com/ekg/vcflib#vcfparsealts), 
[vcfprimers](https://github.com/ekg/vcflib#vcfprimers), 
[vcfqual2info](https://github.com/ekg/vcflib#vcfqual2info), 
[vcfrandom](https://github.com/ekg/vcflib#vcfrandom), 
[vcfrandomsample](https://github.com/ekg/vcflib#vcfrandomsample), 
[vcfremap](https://github.com/ekg/vcflib#vcfremap), 
[vcfremoveaberrantgenotypes](https://github.com/ekg/vcflib#vcfremoveaberrantgenotypes), 
[vcfremovesamples](https://github.com/ekg/vcflib#vcfremovesamples), 
[vcfroc](https://github.com/ekg/vcflib#vcfroc), 
[vcfsample2info](https://github.com/ekg/vcflib#vcfsample2info), 
[vcfsamplediff](https://github.com/ekg/vcflib#vcfsamplediff), 
[vcfsamplenames](https://github.com/ekg/vcflib#vcfsamplenames), 
[vcfsitesummarize](https://github.com/ekg/vcflib#vcfsitesummarize), 
[vcfstats](https://github.com/ekg/vcflib#vcfstats), 
[vcfstreamsort](https://github.com/ekg/vcflib#vcfstreamsort), 
[vcfuniq](https://github.com/ekg/vcflib#vcfuniq), 
[vcfuniqalleles](https://github.com/ekg/vcflib#vcfuniqalleles)

There may be multiple versions of vcflib available. An easy way of selecting the
version is to use [modules](/apps/modules.html). To see the modules
available, type



```

module avail vcflib 

```

To select a module use



```

module load vcflib/[version]

```

where `[version]` is the version of choice.


### Environment variables set


* `$PATH`
* `$CPATH`


### Web sites


* [GitHub](https://github.com/ekg/vcflib#vcfroc)




Interactive job
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **vcf=/fdb/GATK\_resource\_bundle/hg19-2.8/CEUTrio.HiSeq.WGS.b37.bestPractices.hg19.vcf.gz**

[user@cn3144 ~]$ **vcfsamplenames $vcf**
NA12878
NA12891
NA12892

[user@cn3144 ~]$ **zcat $vcf | vcfcountalleles**
12935193

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```



Batch job on Biowulf
Create a batch script similar to the following example:



```

#! /bin/bash

function fail() {
    echo "$@" >&2
    exit 1
}

rb=/fdb/GATK_resource_bundle/hg19-2.8
module load vcflib || fail "could not load vcflib module"
module load samtools/1.2 || fail "could not load samtools module"
tabix -h ${rb}/CEUTrio.HiSeq.WGS.b37.bestPractices.hg19.vcf.gz chr1:1-100000 \
 | vcf2tsv > CEUTrio.out

```


Submit to the queue with [sbatch](/docs/userguide.html):

```

b2$ **sbatch vcf2tsv.sh**

```



Swarm of jobs on Biowulf
Create a swarm command file similar to the following example:



```

vcfannotate -b enhancers.bed -k enh sample1.vcf > sample1_anno.vcf
vcfannotate -b enhancers.bed -k enh sample2.vcf > sample2_anno.vcf

```

And submit to the queue with [swarm](/apps/swarm.html)



```

b2$ **swarm -f vcfannotate.swarm -g 5**

```





