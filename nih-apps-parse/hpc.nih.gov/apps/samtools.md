

document.querySelector('title').textContent = 'samtools on Biowulf';
samtools on Biowulf


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


Samtools is a suite of applications for processing high throughput sequencing data:


* *`samtools`* is used for working with SAM, BAM, and CRAM files containing aligned
 sequences. It can also be used to index fasta files.
* *`bcftools`* is used for working with BCF2, VCF, and gVCF files containing
 variant calls.
* *`bgzip`* is a block compression/decompression utility. `bgzip` compressed
 files can be decompressed with `gunzip`. However, they are in fact made up
 of a number of compressed blocks which facilitates indexing and random access.
* *`tabix`* is used to index `bgzip` compressed tabular data. It can create `.tbi` and
 `.csi` format indices.
* *`wgsim`* is used to simulate NGS reads.
* various conversion utilities (`bowtie2sam`, `novo2sam`, `blast2sam`, `psl2sam`, ...)
* *`htslib`* is a library for reading and writing the formats mentioned above. `samtools`
 and `bcftools` are based on `htslib`


Each of the main tools and file formats have their own manpages
(e. g. `man samtools` or `man sam`).



### Notes on releases


For more detailed release notes see the GitHub release pages for
[htslib](https://github.com/samtools/htslib/releases) and
[samtools](https://github.com/samtools/samtools/releases).


* There have been a number of visible changes with the transition to
the > 1.0 versions based on `htslib`.
* Accessing files via https in addition to http and ftp became available in samtools 1.3. Use
URLs like `http://...`, `https://...`, or
`ftp://...`
* Accessing files stored in Amazon S3 became available in samtools
1.3. Use pseudo URLs in the format `s3://bucket/...`
* Accessing files stored in Google cloud storage became available in
samtools 1.4. Use pseudo URLs in the format `gs://bucket/...`
* Accessing files in the HPC object store is possible since version 1.5 with
a local patch.
* Version 1.10 is currently not able to access the HPC object store


### References:


* Heng Li et al. *The Sequence alignment/map (SAM) format and SAMtools.* 
 Bioinformatics 2009, 25:2078-2079. 
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/19505943) | 
 [PMC](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2723002/) | 
 [Journal](http://bioinformatics.oxfordjournals.org/content/25/16/2078.long)
* Heng Li. *A statistical framework for SNP calling, mutation discovery, 
 association mapping and population genetical parameter estimation from 
 sequencing data*. 
 Bioinformatics 2011, 27:2987-2993.
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/21903627) | 
 [PMC](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3198575/) | 
 [Journal](http://bioinformatics.oxfordjournals.org/content/27/21/2987.long)


Documentation
* [Home page](http://www.htslib.org/)
* [Manual](http://www.htslib.org/doc/#manuals)
* [GitHub](https://github.com/samtools/)


Important Notes
* Module Name: samtools (see [the modules page](/apps/modules.html) for more information)
* Samtools is a multithreaded application
* Example data can be found in `$SAMTOOLS_TEST_DATA`


If the memory per thread is set to a low value for example
by not providing a unit [K/M/G], a potentially **very large**
number of temporary files will be created. Please take care to properly set
this value. Starting with version 1.4 there is a minimum amount of
memory per thread. Check the manpage for details.


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=8g --cpus-per-task=4 --gres=lscratch:30**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load samtools**
[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144 ~]$ **cp $SAMTOOLS\_TEST\_DATA/\* .**
[user@cn3144 ~]$ **ls -l**
-rw-r--r-- 1 user group 45M Jan 18 16:15 mm10_500k_pe.bam
[user@cn3144 ~]$ # show the header
[user@cn3144 ~]$ **samtools view -H mm10\_500k\_pe.bam**
@SQ     SN:chr10        LN:130694993
@SQ     SN:chr11        LN:122082543
@SQ     SN:chr12        LN:120129022
@SQ     SN:chr13        LN:120421639
...snip...

[user@cn3144 ~]$ # sort using lscratch for temp files
[user@cn3144 ~]$ **samtools sort -O bam -o mm10\_500k\_pe\_sorted.bam \
 -T /lscratch/$SLURM\_JOB\_ID/mm10\_500k\_pe \
 -@4 -m 1800M mm10\_500k\_pe.bam**
[user@cn3144 ~]$ # index the sorted file
[user@cn3144 ~]$ **samtools index mm10\_500k\_pe\_sorted.bam**

[user@cn3144 ~]$

```

### Accessing data stored on Amazon AWS and Google cloud storage


Access public data via https or on Amazon S3 and Google cloud storage
(GCS). Only version 1.4 and newer support both protocols:



```

[user@cn3144 ~]$ **module load samtools**
[user@cn3144 ~]$ **samtools view \
https://storage.googleapis.com/genomics-public-data/platinum-genomes/bam/NA12877\_S1.bam \
chr20:100000-100001 | head -5**
[M::test_and_fetch] downloading file 'https://storage.googleapis.com/genomics-public-data/platinum-genom es/bam/NA12877_S1.bam.bai' to local directory
...snip...
[user@cn3144 ~]$ **export HTS\_S3\_V2=1**  ### this is only required for 1.10 for the public buckets
[user@cn3144 ~]$ **samtools view \
s3://1000genomes/phase1/data/NA12878/exome\_alignment/NA12878.mapped.illumina.mosaik.CEU.exome.20110411.bam \
20:100000-100001**
[M::test_and_fetch] downloading file 's3://1000genomes/phase1/data/NA12878/exome_alignment/NA12878.mappe d.illumina.mosaik.CEU.exome.20110411.bam.bai' to local directory
...snip...
[user@cn3144 ~]$ **unset GCS\_OAUTH\_TOKEN**
[user@cn3144 ~]$ **samtools view \
gs://genomics-public-data/platinum-genomes/bam/NA12877\_S1.bam \
chr20:100000-100001 | head -5**
...snip...

```

Note that if you have the environment variable
`GCS_OAUTH_TOKEN` set the last command will fail with a
permissions error. This variable is, however, needed to access data in
your private GCS buckets.


Now, upload data to your private GCS bucket and access it with
samtools. This assumes that you have GCS set up and the google-cloud-sdk
set up and initiated. In the example below, the bucket is named 
`test-2050be34`. One complication is that the python
environments required for Google cloud SDK also contain a samtools
version which may be older than 1.4 and may not work. If you get
'Protocol not supported' errors check that a samtools >= 1.4 is 
first on your path.



```

[user@cn3144 ~]$ **module load google-cloud-sdk**
[user@cn3144 ~]$ # upload some data to your bucket
[user@cn3144 ~]$ **gsutil cp gcat\_set\_053.bam\* gs://test-2050be34/**
[user@cn3144 ~]$ # create a application token with gcloud. This should only be
      # necessary once; this involves a sign-in on the web
[user@cn3144 ~]$ **gcloud auth application-default login**
[user@cn3144 ~]$ # save token to environment variable for samtools.
[user@cn3144 ~]$ **export GCS\_OAUTH\_TOKEN=$(gcloud auth application-default print-access-token)**
[user@cn3144 ~]$ # make sure the python/bin/samtools does not interfere
[user@cn3144 ~]$ **module unload google-cloud-sdk**
[user@cn3144 ~]$ **samtools view gs://test-2050be34/gcat\_set\_053.bam chr20:100000-200001 | head -5**

```

For S3, samtools will try to obtain tokens from either the URL,
environment variables, or ~/.aws/credentials or ~/.s3cfg to access
private buckets. Profiles can be used in which case URLs take the format
`s3://profile@bucket/...`. 

### Accessing data on the HPC object store


Currently not available for 1.10


Samtools starting with version 1.5 can now access data in our local
object store with some preparation. First, let's create or modify a 
`~/.s3cfg` file to be able to access our object store.
Note that the settings are kept in a named profile (obj) rather than the default 
profile. This helps to avoid interfering with the commands above. However,
you may choose to make the local object store the default, in which case the
profile name does not have to be included in the s3 URLs. The important
settings are `access_key, secret_key, host_base, and url_mode`.



```

node$ **cat >> $HOME/.s3cfg <<EOF**
[obj]
access_key = [...your key...]
secret_key = [...your secret...]
host_base = os2access1
host_bucket = os2access1/$USER
url_mode = path
proxy_host = dtn03-e0
proxy_port = 3128
use_https = False

encrypt = False
server_side_encryption = False
human_readable_sizes = True
website_endpoint = http://os2access1/$USER
**EOF**

```

Now copy a bam file and it's index to the object store



```

[user@cn3144 ~]$ **obj\_put gcat\_set\_053.bam**
[user@cn3144 ~]$ **obj\_put gcat\_set\_053.bam.bai**

```

And access them directly with samtools. Note that the `obj@` is
only required if your HPC object store settings are not the default in
`~/.s3cfg`. Note also that for now the `+http`
is required to force samtools to use http. This may change in the future.



```

[user@cn3144 ~]$ **samtools view -H s3+http://obj@user/gcat\_set\_053.bam**
@HD     VN:1.3  SO:coordinate
@SQ     SN:chrM LN:16571
@SQ     SN:chr1 LN:249250621
@SQ     SN:chr2 LN:243199373
@SQ     SN:chr3 LN:198022430
[...snip...]
[user@cn3144 ~]$ **samtools view s3+http://obj@user/gcat\_set\_053.bam chr1:2000000-2000100**
11V6WR1:111:D1375ACXX:1:2304:6028:48213 99      chr1    1999930 60      100M    = 
[...snip...]

```

End the interactive session



```

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. samtools.sh), which uses the input file 'samtools.in'. For example:



```

#!/bin/bash
set -e

module load samtools 
module load gnuplot 
cd /data/$USER/test_data
samtools stats bam/gcat_set_025.bam > bam/gcat_set_025.bam.stats
plot-bamstats -p  bam/gcat_set_025_stats/ bam/gcat_set_025.bam.stats

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] samtools.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. samtools.swarm). For example:



```

samtools rmdup bam/read1_250k_sorted.bam bam/read1_250k_sorted_rmdup.bam
samtools rmdup bam/read1_500k_sorted.bam bam/read1_500k_sorted_rmdup.bam

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f samtools.swarm [-g #] [-t #] --module samtools
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module samtools  Loads the samtools module for each subjob in the swarm 
 | |
 | |
 | |








