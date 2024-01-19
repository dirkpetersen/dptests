

document.querySelector('title').textContent = ' SRA-Toolkit';
 SRA-Toolkit


|  |
| --- |
| 
Quick Links
[SRA Source Repositories](#source)
[Estimating space requirements](#space)
[Downloading from NCBI](#ncbi)
[Batch job on Biowulf](#batch)
[Using Swarm](#swarm)
[Interactive jobs](#int)
[dbGaP downloads](#dbgap)
[Executables](#exec)
[Documentation](#doc)
[Configuring SRA-Toolkit](#config)
[Troubleshooting](#trouble)
 |


The [NCBI SRA SDK](http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=std) generates loading and dumping tools with their
respective libraries for building new and accessing existing runs. 



Documentation
* [NCBI SRA Download Guide](http://www.ncbi.nlm.nih.gov/books/NBK242621/)* [SRA Toolkit documentation](http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc)
* [SRA File Formats Guide](ftp://ftp.era.ebi.ac.uk/meta/doc/sra_1_1/SRA_File_Formats_Guide.pdf)
* **Command line help:** Type the command followed by '-h'
* [fasterq-dump guide](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump)


Important Notes
* Module Name: sratoolkit (see [the modules page](/apps/modules.html) for more information)
* fastq-dump is being deprecated. Use fasterq-dump instead -- it is *much* faster and more efficient.
* fasterq-dump uses temporary space while downloading, so you must make sure you have enough space
* Do not run more than the default 6 threads on Helix.
* To run trimgalore/cutadapt/trinity on these files, the quality header needs to be changed, e.g.

```

sed -r 's/(^[\@\+]SRR\S+)/\1\/1/' SRR10724344_1.filter.fastq
sed -r 's/(^[\@\+]SRR\S+)/\1\/2/' SRR10724344_2.filter.fastq 

```
* fasterq-dump requires tmp space during the download. This temporary directory will use approximately the size of the final output file. On Biowulf, the SRAtoolkit module is set up
to use local disk as the temporary directory. Therefore, if running SRAtoolkit on Biowulf, you must allocate local disk as in the examples below.




SRA Source Repositories

 SRA Data currently reside in 3 NIH repositories:

* NCBI - Bethesda and Sterling
* Amazon Web Services (= 'Amazon cloud' = AWS)
 * Google Cloud Platform (GCP)



Two versions of the data exist: the original (raw) submission, and a normalized (extract, transform, load [ETL]) version. NCBI maintains only ETL data online, while AWS and GCP have both ETL and original submission format. Users who want access to the original bams can only get them from AWS or GCP today.

In the case of ETL data, Sratoolkit tools on Biowulf will always pull from NCBI, because it is obviously nearer and there are no fees.
Most sratoolkit tools such as fasterq-dump will pull ETL data from NCBI. 

NCBI is moving the source bam files to cold cloud storage, so to reliably pull down the source files, users should use the Cloud Data Deliver Service (CDDS) as described at <https://ncbiinsights.ncbi.nlm.nih.gov/2021/09/23/sra-cloud-bucket/>

 If requesting "original submission" files in bam or cram or some other format, they can ONLY be obtained from AWS or GCP and will require that the user provide a cloud-billing account to pay for egress charges. See <https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration> and <https://github.com/ncbi/sra-tools/wiki/04.-Cloud-Credentials>. The user needs to establish account information, register it with the toolkit, and authorize the toolkit to pass this information to AWS or GCP to pay for egress charges.

If you attempt to download non-ETL SRA data from AWS or GCP without the account information, you will see an error message along these lines:

```

Bucket is requester pays bucket but no user project provided.

```


Errors during downloads

It is not unusual for users to get errors while downloading SRA data with prefetch, fasterq-dump, or hisat2, because many people are constantly downloading data and the servers can get overwhelmed. Please see the NCBI SRA page 
[Connection Timeouts](https://github.com/ncbi/sra-tools/wiki/06.-Connection-Timeouts)

Estimating space requirements

fasterq-dump takes significantly more space than the old fastq-dump, as it requires temporary space in addition to the final output. As a rule of thumb, the [fasterq-dump guide](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump) suggests getting the size of the accession using 'vdb-dump', then estimating 7x for the output and 6x for the temp files. For example: 

```

helix% **vdb-dump --info SRR2048331**
acc    : SRR2048331
path   : https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR2048331/SRR2048331.2
**size : 657,343,309**
type   : Table
platf  : SRA_PLATFORM_ILLUMINA
SEQ    : 16,600,251
SCHEMA : NCBI:SRA:Illumina:tbl:q1:v2#1.1
TIME   : 0x0000000056644e79 (12/06/2015 10:04)
FMT    : Fastq
FMTVER : 2.5.4
LDR    : fastq-load.2.5.4
LDRVER : 2.5.4
LDRDATE: Sep 16 2015 (9/16/2015 0:0)

```

Based on the third line, you should have 650 MB \* 7 = 4550 MB =~ 4.5 GB for the tmp files, and 4 GB for the output file(s). It is also recommended that the output file and temporary files be on different filesystems, as in the examples below. 


Downloading data from SRA

You can download SRA fastq files using the fasterq-dump tool, which will download the fastq file into your current working directory by default. 
(Note: the old fastq-dump is being deprecated). During the download, a temporary directory will be created in the location specified by the -t flag (in the example below, in /scratch/$USER) 
that will get deleted after the download is complete.

For example, on Helix, the interactive data transfer system, you can download as in the example below. To download on Biowulf, **don't run on the Biowulf login node**; use a batch job or interactive job instead.



```

[USER@helix]$ **mkdir /data/$USER/sra**     

[USER@helix]$  **module load sratoolkit**

# Note: don't download to /data/$USER, use a subdirectory like /data/$USER/sra instead
[USER@helix]$ **fasterq-dump -p -t /scratch/$USER -O /data/$USER/sra SRR2048331**    
join   :|-------------------------------------------------- 100.00%
concat :|-------------------------------------------------- 100.00%
spots read      : 16,600,251
reads read      : 16,600,251
reads written   : 16,600,251

```

/scratch is not accessible from the Biowulf compute nodes. On a Biowulf interactive session, you should allocate local disk and use that instead of /scratch as in the [example below](#int).

Submitting a single batch job
1. Create a script file similar to the one below.

 

```
#!/bin/bash 

mkdir -p  /data/$USER/sra
module load sratoolkit
fasterq-dump  -t /lscratch/$SLURM_JOBID  -O /data/$USER/sra SRR2048331
sam-dump SRR2048331 > SRR2048331.sam
....
....
```

2. Submit the script on biowulf: 
 
 
```
[biowulf]$ sbatch --gres=lscratch:30  --cpus-per-task=6  myscript
```


Note: this job allocates 30 GB of local disk (--gres=lscratch:30) and then uses the flag -t /lscratch/$SLURM\_JOBID to write temporary files to local disk. 
If you do not allocate local disk and use the -t flag, the temporary files will be written to the current working directory. It is more efficient for your job and for the
system as a whole if you use local disk. See below: 



|  |  |  |  |  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **Command  **TMPDIR  **Output Directory  **Time
| time fasterq-dump -t /lscratch/$SLURM\_JOBID SRR2048331  local disk on Biowulf node  /data/$USER/sra  49 seconds
| time fasterq-dump SRR2048331  /data/$USER/sra  /data/$USER/sra  68 seconds

 | | | |
 | | | |** |** |** |** |



Using Swarm

**NOTE**: The SRA Toolkit executables use random access to read input files. Because of this, users with data located on
GPFS filesystems will see significant slowdowns in their jobs. For SRA data (including dbGaP data) it is best to first copy the input files to a local
/lscratch/$SLURM\_JOBID directory, work on the data in that directory, and copy the results back at the end of the job, as in the example below. See
the section on [using local disk](https://hpc.nih.gov/docs/userguide.html#local) in the Biowulf User Guide.





Using the 'swarm' utility, one can submit many jobs to the cluster to run concurrently. 


Set up a swarm command file (eg /data/username/cmdfile). Here is a sample file that downloads SRA data using fasterq-dump



```
# run fasterq-dump to download the data, then process further, then copy results back to /data

fasterq-dump --aligned --table PRIMARY_ALIGNMENT -O /lscratch/$SLURM_JOBID SRR1234 ; some_command ; \
  cp -R /lscratch/$SLURM_JOBID/some_files  /data/$USER/myoutputdir/
fasterq-dump --aligned --table PRIMARY_ALIGNMENT -O /lscratch/$SLURM_JOBID SRR3456 ; some_command; \
  cp -R /lscratch/$SLURM_JOBID/some_files  /data/$USER/myoutputdir/

```


If you have previously downloaded SRA data into your own directory, you can copy those files to local scratch on the node, process them there, then copy the output back to your /data area. Sample swarm command file:

```

# copy files from /data, run fasterq-dump and other commands, then copy output back to /data
cp /data/user/path/to/SRR1234.sra /lscratch/$SLURM_JOBID; \
  fasterq-dump --aligned --table PRIMARY_ALIGNMENT -O /lscratch/$SLURM_JOBID /lscratch/$SLURM_JOBID/SRR1234.sra ; \
  some_other_command ; \
  cp -R /lscratch/$SLURM_JOBID/some_files /data/$USER/myoutputdir/
cp /data/user/path/to/SRR56789.sra /lscratch/$SLURM_JOBID; \
  fasterq-dump --aligned --table PRIMARY_ALIGNMENT -O /lscratch/$SLURM_JOBID /lscratch/$SLURM_JOBID/SRR56789.sra ; \
  some_other_command ; \
  cp -R /lscratch/$SLURM_JOBID/some_files /data/$USER/myoutputdir/

[....]
```

The **--gres=lscratch:*N*** must be included
in the swarm commands to allocate local disk on the node. 
For example, to allocate 100GB of scratch space and 4GB of memory:



```
$ swarm -f cmdfile --module sratoolkit --gres=lscratch:100 -g 4  -t 6
```

For more information regarding running swarm, see <swarm.html>


Running an interactive job

Allocate an interactive session and run the interactive job there.



```
[biowulf]$ **sinteractive --gres=lscratch:20 --cpus-per-task=6**
salloc.exe: Granted job allocation 789523
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0135 are ready for job

[cn0135]$ **mkdir -p /data/$USER/sra**

[cn0135]$ **module load sratoolkit**

[cn0135]$ **fasterq-dump -t /lscratch/$SLURM\_JOBID SRR2048331 -O /data/$USER/sra**  

[cn0135]$ **exit**
salloc.exe: Job allocation 789523 has been revoked.
[biowulf]$
```





dbGaP download

NCBI's database of Genotypes and Phenotypes (dbGaP) was developed to archive and distribute the data and results from studies 
that have investigated the interaction of genotype and phenotype in Humans. Most dbGaP data is controlled-access. 
[Documentation for downloading dbGap data](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=dbgap_use).

**IMPORTANT:** as of v3.0.0, dbGaP downloads do not work via the proxy servers, and therefore can only
be run on Helix which has a direct connection to the internet. (This will be fixed in the next release) This has two implications: 
* Helix is a single shared system, and each user is limited to 4-6 CPUs at any time. Therefore you cannot run a swarm of download commands on Helix. 
If you have several dbGap accessions to download, they should be done serially as in the example script below.
* There is no /lscratch on Helix. Instead you can use /scratch.



If you are having problems with dbGaP downloads, please try this test download. It accesses a copy of public 1000 Genomes data at NCBI. This is to confirm whether it is a general problem, or specific to your configuration, or specific to the accessions you are trying to download.


```

helix% **module load sratoolkit**

helix% **prefetch --ngc /usr/local/apps/sratoolkit/prj\_phs710EA\_test.ngc \
 -O /data/$USER/test-dbgap SRR1219902**

```



You should see two files called SRR1219902\_dbGaP-0.sra and SRR1219902\_dbGaP-0.sra.vdbcache appear in /data/$USER/test-dbgap/

To download several dbGap accessions on Helix, set up a script along the following lines:

```

#!/bin/bash

module load sratoolkit
prefetch   --ngc my_ngc_file   -O /data/$USER/mydir  SRR111111
prefetch   --ngc my_ngc_file   -O /data/$USER/mydir  SRR111112
etc.

```

or 

```

#!/bin/bash

module load sratoolkit
fasterq-dump --ngc my_ngc_file -t /scratch/$USER -O /data/$USER/mydir SRR111111
fasterq-dump --ngc my_ngc_file -t /scratch/$USER -O /data/$USER/mydir SRR111112
etc

```

Run this script on Helix via

```

helix% **bash my\_download\_script**

```


Executables

As of v 3.0.0, the SRAToolkit contains the following executables:

```

abi-dump          fastq-dump       prefetch       sra-search     vdb-config
abi-load          fastq-load       rcexplain      sra-sort       vdb-copy
align-info        helicos-load     remote-fuser   sra-sort-cg    vdb-decrypt
bam-load          illumina-dump    sam-dump       sra-stat       vdb-dump
blastn_vdb        illumina-load    sff-dump                      vdb-encrypt
cache-mgr         kar              sff-load       sratools       vdb-lock
cg-load           kdbmeta          sra-blastn     srf-load       vdb-passwd
dump-ref-fasta    latf-load        srapath                       vdb-unlock
fasterq-dump      pacbio-load      sra-pileup     test-sra       vdb-validate

```


Configuring SRA-Toolkit On Helix/Biowulf

By default, the SRA Toolkit installed on Biowulf is set up to use the central Biowulf configuration file, which is set up to NOT maintain a local cache of SRA data. After discussion with NCBI SRA
developers, it was decided that this was the most appropriate setup for most users on Biowulf. The [hisat](/apps/hisat.html) program can automatically download
SRA data as needed.

In some cases, users may want to download SRA data and retain a copy. To download using NCBI's 'prefetch' tool, you would need to set up your own configuration file for the NCBI SRA toolkit. Use the command 
vdb-config to set up a directory for downloading. In the following example, the vdb-config utility is used to set up /data/$USER/sra-data as the local
repository for downloading SRA data. Remember that /home/$USER is limited to a quota of 16 GB, so it is best to direct your downloaded SRA data to /data/$USER. 

Sample session: user input in bold.


```
[user@biowulf ~]$ **vdb-config --interactive --interactive-mode textual**
     vdb-config interactive

  data source

   NCBI SRA: enabled (recommended) (1)

   site    : enabled (recommended) (2)


  local workspaces

  Open Access Data
not cached (not recommended) (3)
location: '' (4)


To cancel and exit      : Press 
To update and continue : Enter corresponding symbol and Press 

Your choice > **4**

Path to Public Repository:

Enter the new path and Press 
Press  to accept the path
> **/data/user/sra-data**

Changing root path to 'Public' repository to '/data/user/sra-data'

 vdb-config interactive

 data source

 NCBI SRA: enabled (recommended) (1)

 site : enabled (recommended) (2)


 local workspaces

 Open Access Data
not cached (not recommended) (3)
location: '/data/user/sra-data' (4)


To cancel and exit : Press 
To save changes and exit: Enter Y and Press 
To update and continue : Enter corresponding symbol and Press 

Your choice > **3**

Enabling user repository caching...

 vdb-config interactive

 data source

 NCBI SRA: enabled (recommended) (1)

 site : enabled (recommended) (2)


 local workspaces

 Open Access Data
cached (recommended) (3)
location: '/data/user/sra-data' (4)

To cancel and exit : Press 
To save changes and exit: Enter Y and Press 
To update and continue : Enter corresponding symbol and Press 

Your choice > **y**
Saving...
Exiting...

[biowulf]$ 
```

For more information about encrypted data, please see
[Protected Data Usage Guide](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=dbgap_use) at NCBI.



Troubleshooting

* For dbGaP downloads, try the test download on Helix first: 

```

helix% module load sratoolkit
helix% prefetch --ngc /usr/local/apps/sratoolkit/prj_phs710EA_test.ngc -O /data/$USER SRR1219902

```

If you get an error, try renaming your ~/.ncbi configuration directory.

```

helix% mv ~/.ncbi ~/.ncbi.OLD

```

If you still get an error downloading dbGaP data on helix, let us know by sending email to staff@hpc.nih.gov.

* If you are getting this error while using fastq-dump and fasterq-dump: 

```

Failed to call external services

```

try this:

```

% **mv ~/.ncbi ~/.ncbi.OLD**

```

This will move your existing SRA toolkit configuration out of the way, to test whether some setting in your configuration is causing the problem.

* If you get the error

```

fasterq-dump was killed (signal 9 SIGKILL)

```

you are likely running fasterq-dump on a compute node and have not allocated enough memory. Try increasing the memory
allocation by 10-20 GB and rerun. Or run your download on Helix.










































































