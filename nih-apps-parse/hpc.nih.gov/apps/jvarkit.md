

document.querySelector('title').textContent = 'jvarkit on Biowulf';
jvarkit on Biowulf


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



 Java tools for bioinformatics. Tools are not described in detail here see
 the author's [documentation](http://lindenb.github.io/jvarkit/).
 Note that all tools are compiled though some of them are pretty lab specific
 or deprecated.


### References:


* Lindenbaum, Pierre. *JVarkit: java-based utilities for Bioinformatics.*. 
 <http://dx.doi.org/10.6084/m9.figshare.1425030>



Documentation
* jvarkit on [GitHub](https://github.com/lindenb/jvarkit)


Important Notes
* Module Name: jvarkit (see [the modules page](/apps/modules.html) 
 for more information)
* This is a collection of jar files all located in `$JVARKIT_JARPATH`
* Example files in `$JVARKIT_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=6g --gres=lscratch:20**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load jvarkit**
[user@cn3144]$ **cp -L ${JVARKIT\_TEST\_DATA:-none}/\* .**
[user@cn3144]$ **java -jar $JVARKIT\_JARPATH/backlocate.jar --help**
Usage: backlocate [options] Files
  Options:
  * -g, --gtf
      A GTF (General Transfer Format) file. See
      https://www.ensembl.org/info/website/upload/gff.html . Please note that
      CDS are only detected if a start and stop codons are defined.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
    -p, --printSeq
      print mRNA & protein sequences
      Default: false
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools
      faidx and with picard CreateSequenceDictionary
    --version
      print version and exit


[user@cn3144]$ **gc=/fdb/GENCODE/Gencode\_human/release\_35**
[user@cn3144]$ **gtf=$gc/gencode.v35.primary\_assembly.annotation.gtf** 
[user@cn3144]$ **genome=$gc/GRCh38.primary\_assembly.genome.fa**

[user@cn3144]$ **echo -e "NOTCH2\tPro1090M\tInteresting" \
 | java -jar $JVARKIT\_JARPATH/backlocate.jar \
 --gtf $gtf -R $genome \
 | grep -v "##" \
 | java -jar $JVARKIT\_JARPATH/prettytable.jar** 
+------------+-----+--------------+-----+-----------------+-------------------+-------------------+---------------+---------------+------------+----------------------+-------------+------------+-------------------+--------------------------+----------+-----------------+
| #User.Gene | AA1 | petide.pos.1 | AA2 | transcript.name | transcript.id     | transcript.strand | transcript.AA | index0.in.rna | wild.codon | potential.var.codons | base.in.rna | chromosome | index0.in.genomic | exon                     | messages | extra.user.data |
+------------+-----+--------------+-----+-----------------+-------------------+-------------------+---------------+---------------+------------+----------------------+-------------+------------+-------------------+--------------------------+----------+-----------------+
| NOTCH2     | Pro | 1090         | Met | NOTCH2          | ENST00000256646.7 | -                 | P             | 3267          | CCA        | .                    | C           | chr1       | 119937925         | ENST00000256646.7.Exon20 | .        | Interesting     |
| NOTCH2     | Pro | 1090         | Met | NOTCH2          | ENST00000256646.7 | -                 | P             | 3268          | CCA        | .                    | C           | chr1       | 119937924         | ENST00000256646.7.Exon20 | .        | Interesting     |
| NOTCH2     | Pro | 1090         | Met | NOTCH2          | ENST00000256646.7 | -                 | P             | 3269          | CCA        | .                    | A           | chr1       | 119937923         | ENST00000256646.7.Exon20 | .        | Interesting     |
+------------+-----+--------------+-----+-----------------+-------------------+-------------------+---------------+---------------+------------+----------------------+-------------+------------+-------------------+--------------------------+----------+-----------------+

[user@cn3144]$ **paste <(gunzip -c r1.fastq.gz | paste - - - -) <(gunzip -c r2.fastq.gz | paste - - - -)\
 | tr "\t" "\n" \
 | java -jar $JVARKIT\_JARPATH/fastqjs.jar \
 -i -e 'pair.get(0).getReadString().contains("GTTTAAAC") && pair.get(1).getReadString().contains("GTTTAAAC") ' \ 
 > example.fastq** 
[user@cn3144]$ **cat <<\_\_EOF\_\_ > filter.js**
// prints a VARIATION if two samples at least have a DP<200 
function myfilterFunction()
    {
    var samples=header.genotypeSamples;
    var countOkDp=0;


    for(var i=0; i< samples.size();++i)
        {
        var sampleName=samples.get(i);
        if(! variant.hasGenotype(sampleName)) continue;
        var genotype = variant.genotypes.get(sampleName);
        if( ! genotype.hasDP()) continue;
        var dp= genotype.getDP();
        if(dp < 200 ) countOkDp++;
        }
    return (countOkDp>2)
    }
myfilterFunction();
**\_\_EOF\_\_**   
[user@cn3144]$ **java -jar $JVARKIT\_JARPATH/vcffilterjs.jar -f filter.js gatk.vcf**
...
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  BLANK   NA12878 NA12891 NA12892 NA19238  NA19239 NA19240
chr22   42526449        .       T       A       151.47  .       AC=1;AF=0.071;AN=14;BaseQRankSum=2.662;DP=1226;DS;Dels=0.00;FS=0.000;HRun=0;HaplotypeScore=41.2083;MQ=240.47;MQ0=0;MQRankSum=0.578;QD=4.89;ReadPosRankSum=3.611   GT:AD:DP:GQ:PL  0/1:23,8:31:99:190,0,694        0/0:188,0:190:99:0,478,5376     0/0:187,0:187:99:0,493,5322      0/0:247,0:249:99:0,634,6728     0/0:185,0:185:99:0,487,5515     0/0:202,0:202:99:0,520,5857      0/0:181,1:182:99:0,440,5362
chr22   42526634        .       T       C       32.60   .       AC=1;AF=0.071;AN=14;BaseQRankSum=1.147;DP=1225;DS;Dels=0.00;FS=0.000;HRun=0;HaplotypeScore=50.0151;MQ=240.65;MQ0=0;MQRankSum=1.151;QD=1.30;ReadPosRankSum=1.276   GT:AD:DP:GQ:PL  0/1:21,4:25:71:71,0,702 0/0:187,2:189:99:0,481,6080     0/0:233,0:233:99:0,667,7351      0/0:230,0:230:99:0,667,7394     0/0:174,1:175:99:0,446,5469     0/0:194,2:196:99:0,498,6239      0/0:174,0:175:99:0,511,5894
chr22   42527793        rs1080989       C       T       3454.66 .       AC=2;AF=0.167;AN=12;BaseQRankSum=-3.007;DB;DP=1074;DS;Dels=0.01;FS=0.000;HRun=1;HaplotypeScore=75.7865;MQ=209.00;MQ0=0;MQRankSum=3.014;QD=9.36;ReadPosRankSum=0.618       GT:AD:DP:GQ:PL  ./.     0/1:72,90:162:99:1699,0,1767    0/1:103,96:202:99:1756,0,2532    0/0:188,0:188:99:0,526,5889     0/0:160,0:160:99:0,457,4983     0/0:197,0:198:99:0,544,6100      0/0:156,0:156:99:0,439,5041

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. jvarkit.sh), which uses the input file 'jvarkit.in'. For example:



```

#!/bin/bash
module load jvarkit/20200713
paste <(gunzip -c r1.fastq.gz | paste - - - -) <(gunzip -c r2.fastq.gz | paste - - - -) \
    | tr "\t" "\n" \
    | java -jar $JVARKIT_JARPATH/fastqjs.jar \
        -i -e 'pair.get(0).getReadString().contains("GTTTAAAC") && pair.get(1).getReadString().contains("GTTTAAAC") ' \ 
    > example.fastq

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --mem=2g jvarkit.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. jvarkit.swarm). For example:



```

echo -e "NOTCH2\tP1090M\tInteresting" | java -jar $JVARKIT_JARPATH/backlocate.jar --gtf hg19.gtf -R hg19.fa > notch_P1090M
echo -e "NOTCH2\tP1090I\tInteresting" | java -jar $JVARKIT_JARPATH/backlocate.jar --gtf hg19.gtf -R hg19.fa > notch_P1090I

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f jvarkit.swarm -g 2 --module jvarkit/20200713
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module jvarkit  Loads the jvarkit module for each subjob in the swarm 
 | |
 | |
 | |








