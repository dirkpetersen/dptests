

document.querySelector('title').textContent = 'fragpipe on Biowulf';
fragpipe on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



FragPipe is a Java Graphical User Interface (GUI) and CLI workflow tool for a suite of computational tools enabling comprehensive analysis of mass spectrometry-based proteomics data. It is powered by MSFragger - an ultrafast proteomic search engine suitable for both conventional and “open” (wide precursor mass tolerance) peptide identification. 



Documentation
* [fragpipe Main Site](https://fragpipe.nesvilab.org/)


Important Notes
GUI use of this application requires a [graphical connection using NX](https://hpc.nih.gov/docs/nx.html)


* Beginning with version 18.0, fragpipe may be run in headless mode. See the [Batch Job](#sbatch) section for more details.
 * Module Name: fragpipe (see [the modules page](/apps/modules.html) for more information)
 * Environment variables set 
	+ FRAGPIPE\_HOME
	+ FRAGPIPE\_TAR* fragpipe saves configuration and session information in the same directory where it is installed. If you want your data to persist across sessions, you need to copy this program to your own space. See below.



Interactive job
The fragpipe GUI must be run in an [Interactive job](/docs/userguide.html#int). The underlying programs that it uses like [msfragger](msfragger.html) and [philosopher](philosopher.html) can be run as batch or swarm jobs.
First, start a [nomachine desktop session](https://hpc.nih.gov/docs/nx.html) on the Biowulf login node. Within your nx session, allocate an [interactive session](/docs/userguide.html#int) and run fragpipe.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load fragpipe**

[user@cn3144 ~]$ **fragpipe**

```


You should now see the Java Graphical User Interface (GUI).

 ![fragpipe image](/images/fragpipe.PNG)



If you run fragpipe from this central location, you will be unable to save data (like newly created pipelines) and recover them in a new session. You can copy fragpipe to your own space using the following commands or similar:


```

[user@cn3144]$ **module load fragpipe**

[user@cn3144]$ **export INSTALL\_DIR=/data/${USER}/fragpipe**

[user@cn3144]$ **mkdir -p $INSTALL\_DIR && cd $INSTALL\_DIR**

[user@cn3144]$ **tar xvf $FRAGPIPE\_TAR**

```


These steps need only be carried out once. In future sessions fragpipe will already be installed in your space. You can now run it using the full path (where <ver> is the version of fragpipe that you copied).


```

[user@cn3144]$ **/data/${USER}/fragpipe/<ver>/bin/fragpipe**

```


Or you can add the directory to your path to open more easily:


```

[user@cn3144]$ **export PATH=/data/${USER}/fragpipe/<ver>/bin:$PATH**

[user@cn3144]$ **fragpipe**

```


In either case, your saved data will now persist across sessions. 
When running fragpipe from your own space, it's still a good idea to load the fragpipe module first to properly prepare your environment. 

Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
FragPipe versions 18.0+ may now be run in [headless mode](https://fragpipe.nesvilab.org/docs/tutorial_headless.html) to submit as part of an interactive or batch job.


Workflow and manifest files should be created using the GUI either in an interactive session or on your own computer.


 ![fragpipe headless config](/images/fragpipe-headless.jpg)


Create a batch input file (e.g. msfragger.sh). For example:



```

#!/bin/bash
set -e
module load fragpipe/18.0
cd /path/to/data/and/conf

fragpipe --headless \
  --workflow <path to workflow file> \
  --manifest <path to manifest file> \
  --workdir <path to result directory>

```



















