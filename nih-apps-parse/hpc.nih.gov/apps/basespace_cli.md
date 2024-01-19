

document.querySelector('title').textContent = ' BaseSpace CLI on Biowulf &amp; Helix';

BaseSpace CLI on Biowulf & Helix



|  |
| --- |
| 
Quick Links
[How to](#howto)
 |



Description

The Illumina BaseSPace Sequence Hub is a cloud based platform for
analyzing data from Illumina sequencers. It directly integrates with
sequecing machines to monitor runs and stream data to BaseSpace.
Predefined pipelines can be used to analyze the data streamed from the
sequencers or uploaded through another mechanism. 




Storage and compute are provided by AWS.




BaseSpace Sequence Hub can be accessed through its web interface as well
as through the command line interface (CLI) described here.



There may be multiple versions of BaseSpace CLI available. An easy way of selecting the
version is to use [modules](/apps/modules.html). To see the modules
available, type



```

module avail basespace_cli 

```

To select a module use



```

module load basespace_cli/[version]

```

where `[version]` is the version of choice.


### Environment variables set


* $PATH


### Documentation


* [Manual](https://help.basespace.illumina.com/articles/descriptive/basespace-cli)
* [Overview](https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview)




How to
The BaseSpace CLI provides an interactive interface with Illumina's
BaseSpace. Actual computing and storage is done in the cloud. Therefore
this will generally be used interactively. For this example we will use
an interactive session. Before using the CLI for the first time,
it will be necessary to authenticate. This will store credentials
necessary to access a BaseSpace account in `$HOME/.basespace`.
Visit the URL provided by `bs authenticate` to create
the required access token



```

biowulf$ **sinteractive**
salloc.exe: Pending job allocation 21758857
[...snip...]
salloc.exe: Nodes cn2623 are ready for job
cn2623$ **module load basespace\_cli**
[+] Loading basespace_cli 0.8.1
cn2623$ **bs authenticate**
please authenticate here:
https://basespace.illumina.com/oauth/device?code=XXXXX
...
Success!

```

Create a project:



```

cn2623$ **bs list projects**
# no projects yet
cn2623$ **bs create project -n "TestProject"**
cn2623$ **bs list projects**
+------------+--------------+
| project id | project name |
+------------+--------------+
| 55555555   | TestProject  |
+------------+--------------+

```

Upload an illumina generated sample.



```

cn2623$ **bs create biosample -n "MySample" -p "TestProject"**
+-------------+-----------+-----------+--------------------------------------+
|    Name     |    Id     | TotalSize |             DateCreated              |
+-------------+-----------+-----------+--------------------------------------+
| TestProject | 256402146 | 0         | 2021-05-04 13:28:29.018859 +0000 UTC |
+-------------+-----------+-----------+--------------------------------------+
cn2623$ **bs upload sample -p "256402146" TestSample\_S1\_L001\_R1\_001.fastq.gz**
Gathering metadata and validating fastq files...
Uploading ...
        TestSample_S1_L001_R1_001.fastq.gz ..... complete 
Uploaded sample with ID: 38322428
Uploaded by #### ####, using BaseSpaceCLI.SampleUpload/0.8.1 v0.8.1 on biowulf.nih.gov
cn2623$ **bs list biosample --project-name TestProject**
+-----------+-------------+
| sample id | sample name |
+-----------+-------------+
| 38322428  | TestSample  |
+-----------+-------------+

```

Some applications are set up by default. Others have to be imported.
The easiest way to do this is based on a previous run of the application.
Here we use the SRA importer application.



```

cn2623$ **bs list application --filter-term=SRA**
+-----------------------------+----------+---------------+
|            Name             |    Id    | VersionNumber |
+-----------------------------+----------+---------------+
| SRA Submission - DEPRECATED | 147147   | 0.0.3         |
| SRA Submission              | 10596586 | 1.0.0         |
| SRA Import                  | 10741731 | 0.0.6         |
+-----------------------------+----------+---------------+

```

Now, launch the SRA import app to import an (Illumina) run from SRA
directly into BaseSpace



```

cn2623$ **bs launch application -i 10741731 -o "sra-id:SRR292678" TestProject**
SRA Import :  (36554570)

```


Please see the
 [Manual](https://help.basespace.illumina.com/articles/descriptive/basespace-cli)
for more detail.







