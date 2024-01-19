

document.querySelector('title').textContent = 'IGV and IGVTools on Biowulf';
IGV and IGVTools on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
 |


The Integrative Genomics Viewer (IGV) is a high-performance visualization tool for interactive exploration of large, integrated datasets. Features include:


* Support for a wide variety of data types including
+ Genetic variation (copy number, loss of heterozygosity, somatic mutations)
+ Gene/microRNA expression
+ RNAi screens
+ Epigenetic data
+ Genomic annotations
+ Sequence Alignments

* Easy access to genomes and datasets hosted by the Broad Institute


### References:


* James T. Robinson, Helga Thorvaldsdóttir, Wendy Winckler, Mitchell Guttman, Eric S. Lander, Gad Getz, Jill P. Mesirov. [Integrative Genomics Viewer. Nature Biotechnology 29, 24–26 (2011)](http://www.nature.com/nbt/journal/v29/n1/abs/nbt.1754.html)
* Helga Thorvaldsdóttir, James T. Robinson, Jill P. Mesirov. [Integrative Genomics Viewer (IGV): high-performance genomics data visualization and exploration. Briefings in Bioinformatics 14, 178-192 (2013).](https://academic.oup.com/bib/article/14/2/178/208453/Integrative-Genomics-Viewer-IGV-high-performance?searchresult=1)


Documentation
* [IGV Main Site](http://software.broadinstitute.org/software/igv/)


Important Notes
This application requires a [graphical connection using NX](/docs/nx.html)


* Module Name: IGV and IGVTools (see [the modules page](/apps/modules.html) for more information)

* environment variables set 
	+ IGV\_HOME
	+ IGV\_JAR
	+ IGV\_JARPATH
	+ IGVTOOLS\_HOME
	+ IGVTOOLS\_JAR
	+ IGVTOOLS\_JARPATH* Mac users: If mouse right-click fail to work when opening IGV in NoMachine, please 'peel' the top right corner of NoMachine window, select 'input', then check 'grab the mouse input'.



**NOTE:** When running on Biowulf, make sure that web proxying is enabled by the **http\_proxy** environment variable.
Otherwise you may get a java error message about being unable to connect to remote hosts. For 
information about web proxies on the Biowulf cluster, see [here](https://hpc.nih.gov/docs/transfer.html#compute ).


Common pitfall
This is a list of common errors for IGV users on the cluster. 



**Failed to load saved session**
This can be fixed by unusing the relative path:

```

View -> Preferences -> General
Unselect "Use relative paths in session files"

```

The cause of this error is IGV can't locate the files using relative paths. 
This appears to be a problem under only some circumstances.


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=8g --gres=lscratch:5**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load IGV IGVTools**
[user@cn3144 ~]$ **igv -m 8g**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

Here are the commandline options:


* **-b [file]** or **--batch=[file]** Immediately run the supplied batch script after start up.
* **-g [genomeId]** or **--genome=[genomeId]** Specify the genome by ID. You can locate the genomeId as the last entry in each line of the genome registry. NEW since 2.1.28: the -g option can specify a path to a .genome or indexed fasta file.
* **-d [URL]** or **--dataServerURL=[URL]** Specify URL to data registry file.
* **-u [URL]** or **--genomeServerURL=[URL]** Specify URL to genome registry file.
* **-p [port]** or **--port=[port]** Specify port IGV opens for port commands.
* **-o [file]** or **--preferences=[file]** A user preference properties file. This can override any properties IGV supports.


The igv command on Helix is wrapped to include two additional options:


* **-m** or **--memory** memory allocated
* **-t** or **--tmpdir** temporary directory


By default, IGV uses 2gb of memory and uses either **/lscratch/$SLURM\_JOBID/igvtmp** if **--gres=lscratch:*N*** is used in a batch job (where ***N*** is the number GB needed for scratch space, where
 **$USER** is the login of the user running IGV. To allocate 20gb of memory, use:



```
igv --memory 20g
```

Input/output for IGV is menu-driven through the GUI.


IGVTools
--------


IGVTools allows commandline utilities for working with ascii file formats.


* **tile** converts a sorted data input file to a binary tiled data (.tdf) file. Supported input file formats: .wig, .cn, .snp, .igv, .gct
* **count** computes average alignment or feature density for over a specified window size across the genome. Supported input file formats: .sam, .bam, .aligned, .sorted.txt, .bed
* **sort** sorts the input file by start position. Supported input file formats: .cn, .igv, .sam, .aligned, and .bed
* **index** creates an index file for an input ascii alignment or feature file. Supported input file formats: .sam, .aligned, .sorted.txt, and .bed


By default, IGVTools uses 5gb of memory and the same temporary directory method as IGV.
 To allocate 20gb of memory, use:



```
igvtools --memory 20g [ additional options for igvtools ]
```

To see what options are available for the commandline IGVTools, just type **igvtools** at the prompt:



```

Program: igvtools. IGV Version 2.3.98 (141)07/25/2017 12:12 AM

Usage: igvtools [command] [options] [input file/dir] [other arguments]

Command: version print the version number
	 sort    sort an alignment file by start position.
	 index   index an alignment file
	 toTDF    convert an input file (cn, gct, wig) to tiled data format (tdf)
	 count   compute coverage density for an alignment file
	 formatexp  center, scale, and log2 normalize an expression file
	 gui      Start the gui
	 help      display this help message, or help on a specific command
	 See http://www.broadinstitute.org/software/igv/igvtools_commandline for more detailed help

```





