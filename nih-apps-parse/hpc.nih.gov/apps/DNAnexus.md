

document.querySelector('title').textContent = 'DNAnexus on Biowulf';
DNAnexus on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive session](#helix)
[Batch job](#batch)
 |



DNAnexus is a cloud-based commercial solution for next-generation sequence analysis and visualization. It has a command-line interface (CLI) which can be used to log in
to the DNAnexus platform, upload and navigate data, and launch analyses. For users who wish to upload/download data directly from the NIH HPC systems, the CLI has been 
installed on Biowulf



Documentation
* [DNAnexus CLI quick start](https://documentation.dnanexus.com/getting-started/cli-quickstart)
* [dx-toolkit Github](https://github.com/dnanexus/dx-toolkit)


Important Notes
* Module Name: DNAnexus (see [the modules page](/apps/modules.html) for more information)
* Requires a DNAnexus account set up before use. Accounts can be set up at [platform.dnanexus.com](https://platform.dnanexus.com/). Small test jobs are 
free but large-scale analyses will require billing to be set up by each user.



Helix, or Biowulf compute nodes


 **DNAnexus cannot be run on the Biowulf login node.** 
Since DNAnexus CLI commands are either doing data transfer or remote jobs, they can be safely run on Helix (the interactive data transfer system where 
scientific applications are not available), or you can allocate an interactive session on Biowulf by typing 'sinteractive'. The example below shows
an interactive DNAnexus session on Helix, but an interactive session on a Biowulf compute node would be identical.


```

helix% **module load DNAnexus**
helix% **dx login**
Acquiring credentials from https://auth.dnanexus.com
Username [user]:
Password:
Using env variable HTTPS_PROXY=http://dtn04-e0:3128 as proxy

Note: Use dx select --level VIEW or dx select --public to select from projects for
which you only have VIEW permissions.

Available projects (CONTRIBUTE or higher):
0) Test (ADMINISTER)

Pick a numbered choice [0]: **0**
Setting current project to: Test
You are now logged in. Your credentials are stored in /home/user/.dnanexus_config and will expire
in 30 days, 0:00:00. Use dx login --timeout to control the expiration date, or dx
logout to end this session.

helix% **dx select**
Available projects (CONTRIBUTE or higher):
0) Test (ADMINISTER)
Pick a numbered choice [0]: **0**
Setting current project to: Test

helix% **dx upload --wait /data/$USER/fastq\_files/8\_2\_61PKFAAXX.fq.gz**
[===========================================================>] Uploaded 257,945,911 of 257,945,911 bytes (100%) /data/CloudData/sean/fastq_files/8_2_61PKFAAXX.fq.gz
ID                  file-FJzQ1b003PyyyjVbBVq84xKK
Class               file
Project             project-FJzPQV003PyqqVx409PJBgZ0
Folder              /
Name                8_2_61PKFAAXX.fq.gz
State               closed
Visibility          visible
Types               -
Properties          -
Tags                -
Outgoing links      -
Created             Thu Aug 23 14:45:40 2018
Created by          user
Last modified       Thu Aug 23 14:47:04 2018
archivalState       "live"
Media type          application/gzip
Size                246.00 MB

helix% **dx find apps --installed**
x Bowtie2 FASTA Indexer (bowtie2_fasta_indexer), v1.3.0
x BWA FASTA Indexer (bwa_fasta_indexer), v2.0.1
x FastQC Reads Quality Control (fastqc), v2.3.0

helix% **dx run fastqc**
Entering interactive mode for input selection.

Input:   Reads (reads)
Class:   file

Enter file ID or path ( twice for compatible files in current directory, '?' for
more options)
reads: **8\_2\_61PKFAAXX.fq.gz**
Select an optional parameter to set by its # (^D or  to finish):

 [0] Custom contaminants (contaminants\_txt)
 [1] Custom adapters (adapters\_txt)
 [2] Custom limits (limits\_txt)
 [3] Input format (format) [default="auto"]
 [4] Kmer size (kmer\_size) [default=7]
 [5] Disable grouping of bases past 50bp? (nogroup) [default=true]
 [6] Extra command-line options (extra\_options)

Optional param #:

Using input JSON:
{
 "reads": {
 "$dnanexus\_link": {
 "project": "project-FJzPQV003PyqqVx409PJBgZ0",
 "id": "file-FJzQ1b003PyyyjVbBVq84xKK"
 }
 }
}

Confirm running the executable with this input [Y/n]: **y**
Calling app-FBFjGpQ987vZKJXvF1zfPx8X with output destination project-FJzPQV003PyqqVx409PJBgZ0:/

Job ID: job-FJzQ6K003PyY8BbgBV4qxXYF
Watch launched job now? [Y/n] **y**

Job Log
-------
Watching job job-FJzQ6K003PyY8BbgBV4qxXYF. Press Ctrl+C to stop.
\* FastQC Reads Quality Control (fastqc:main) (running) job-FJzQ6K003PyY8BbgBV4qxXYF
 user 2018-08-23 14:55:52 (running for 0:00:39)
2018-08-23 14:56:36 FastQC Reads Quality Control INFO Logging initialized (priority)
2018-08-23 14:56:37 FastQC Reads Quality Control INFO CPU: 12% (4 cores) \* Memory: 661/7478MB \* Storage: 79GB free \* Net: 0?/0?MBps
2018-08-23 14:56:37 FastQC Reads Quality Control STDOUT dxpy/0.259.0 (Linux-4.4.0-98-generic-x86\_64-with-Ubuntu-14.04-trusty)
2018-08-23 14:56:38 FastQC Reads Quality Control INFO Downloading bundled file resources.tar.gz
2018-08-23 14:56:39 FastQC Reads Quality Control STDOUT >>> Unpacking resources.tar.gz to /
2018-08-23 14:56:39 FastQC Reads Quality Control INFO Installing apt packages openjdk-7-jre-headless
2018-08-23 14:56:59 FastQC Reads Quality Control STDOUT dxpy/0.259.0 (Linux-4.4.0-98-generic-x86\_64-with-Ubuntu-14.04-trusty)
2018-08-23 14:56:59 FastQC Reads Quality Control STDOUT /usr/sbin/sshd already running.
2018-08-23 14:56:59 FastQC Reads Quality Control STDOUT bash running (job ID job-FJzQ6K003PyY8BbgBV4qxXYF)
2018-08-23 14:57:00 FastQC Reads Quality Control STDERR ++ mark-section 'streaming input data'
2018-08-23 14:57:00 FastQC Reads Quality Control STDERR ++ mkfifo ./8\_2\_61PKFAAXX.fq.gz
2018-08-23 14:57:00 FastQC Reads Quality Control STDERR ++ cat\_pid=3235
2018-08-23 14:57:00 FastQC Reads Quality Control STDERR ++ dx cat '{"$dnanexus\_link": "file-FJzQ1b003PyyyjVbBVq84xKK"}'
2018-08-23 14:57:00 FastQC Reads Quality Control STDERR ++ mark-section 'running FastQC'
2018-08-23 14:57:00 FastQC Reads Quality Control STDERR ++ opts=
2018-08-23 14:57:00 FastQC Reads Quality Control STDERR ++ '[' '' '!=' '' ']'
2018-08-23 14:57:00 FastQC Reads Quality Control STDERR ++ '[' '' '!=' '' ']'
2018-08-23 14:57:00 FastQC Reads Quality Control STDERR ++ '[' '' '!=' '' ']'
2018-08-23 14:57:00 FastQC Reads Quality Control STDERR ++ '[' auto '!=' auto ']'
2018-08-23 14:57:00 FastQC Reads Quality Control STDERR ++ '[' true == true ']'
2018-08-23 14:57:00 FastQC Reads Quality Control STDERR ++ opts=' --nogroup'
2018-08-23 14:57:00 FastQC Reads Quality Control STDERR ++ mkdir results
2018-08-23 14:57:00 FastQC Reads Quality Control STDERR +++ nproc
2018-08-23 14:57:00 FastQC Reads Quality Control STDERR ++ /FastQC/fastqc -t 4 --extract -k 7 -o results --nogroup ./8\_2\_61PKFAAXX.fq.gz
2018-08-23 14:57:01 FastQC Reads Quality Control STDERR Started analysis of 8\_2\_61PKFAAXX.fq.gz
2018-08-23 14:57:30 FastQC Reads Quality Control STDOUT Analysis complete for 8\_2\_61PKFAAXX.fq.gz
2018-08-23 14:57:34 FastQC Reads Quality Control STDERR ++ wait 3235
2018-08-23 14:57:34 FastQC Reads Quality Control STDERR ++ mark-section 'uploading results'
2018-08-23 14:57:34 FastQC Reads Quality Control STDERR ++ mkdir -p /home/dnanexus/out/report\_html/ /home/dnanexus/out/stats\_txt/
2018-08-23 14:57:34 FastQC Reads Quality Control STDERR ++ mv results/8\_2\_61PKFAAXX.fq\_fastqc/fastqc\_report.html /home/dnanexus/out/report\_html/8\_2\_61PKFAAXX.stats-fastqc.html
2018-08-23 14:57:34 FastQC Reads Quality Control STDERR ++ mv results/8\_2\_61PKFAAXX.fq\_fastqc/fastqc\_data.txt /home/dnanexus/out/stats\_txt/8\_2\_61PKFAAXX.stats-fastqc.txt
2018-08-23 14:57:34 FastQC Reads Quality Control STDERR ++ dx-upload-all-outputs
2018-08-23 14:57:34 FastQC Reads Quality Control STDOUT uploading file: /home/dnanexus/out/report\_html/8\_2\_61PKFAAXX.stats-fastqc.html -> /8\_2\_61PKFAAXX.stats-fastqc.html
2018-08-23 14:57:35 FastQC Reads Quality Control STDOUT uploading file: /home/dnanexus/out/stats\_txt/8\_2\_61PKFAAXX.stats-fastqc.txt -> /8\_2\_61PKFAAXX.stats-fastqc.txt
2018-08-23 14:57:35 FastQC Reads Quality Control STDERR ++ mark-success
2018-08-23 14:57:35 FastQC Reads Quality Control STDERR ++ type -t main
2018-08-23 14:57:35 FastQC Reads Quality Control STDERR + [[ '' == \f\u\n\c\t\i\o\n ]]
2018-08-23 14:57:35 FastQC Reads Quality Control STDERR + echo '/home/dnanexus/job-FJzQ6K003PyY8BbgBV4qxXYF: Global scope execution complete. Not invoking entry point function main because it was not found'
2018-08-23 14:57:35 FastQC Reads Quality Control STDERR /home/dnanexus/job-FJzQ6K003PyY8BbgBV4qxXYF: Global scope execution complete. Not invoking entry point function main because it was not found
\* FastQC Reads Quality Control (fastqc:main) (done) job-FJzQ6K003PyY8BbgBV4qxXYF
 user 2018-08-23 14:55:52 (runtime 0:01:07)
 Output: report\_html = file-FJzQ77Q0Gjvj15Vz10y6JGXz
 stats\_txt = file-FJzQ77j0Gjvq6kg010ZzZ6Z0
 
helix% **dx ls**
8\_2\_61PKFAAXX.fq.gz
8\_2\_61PKFAAXX.stats-fastqc.html
8\_2\_61PKFAAXX.stats-fastqc.txt

helix%

```


Batch Job on Biowulf

The initial DNAnexus login cannot be performed as part of a Biowulf batch job. You will need to run 'dx login' on interactively on Helix or a Biowulf compute node. Once logged in, the login will be valid for 30 days. 
During that period it is possible to submit a batch job on Biowulf that runs DNAnexus commands.

Sample batch script

```

#!/bin/bash

module load DNAnexus

indir=/data/$USER/fastq_files
infile=8_2_61PKFAAXX
dx upload --wait ${indir}/${infile}.fq.gz
dx ls
dx run --wait fastqc -ireads=$infile
dx ls
dx download ${infile}.stats-fastqc.txt

```


Submit with:

```

sbatch myjob.sh

```

The status of the job and any errors can be seen in the slurm-\*.out file. DNAnexus will also send an
email as below: 

![](/images/DNAnexus_email.png)




















