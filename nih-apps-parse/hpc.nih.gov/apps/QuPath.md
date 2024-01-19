

document.querySelector('title').textContent = 'QuPath on Biowulf';
QuPath on Biowulf


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



"an open, powerful, flexible, extensible software platform for bioimage analysis"



### References:


* [Bankhead, Peter, et al. "QuPath: Open source software for digital pathology image analysis." *Scientific reports* 7.1 (2017): 1-7.](https://www.nature.com/articles/s41598-017-17204-5)


Documentation
* [QuPath Main Site](https://qupath.github.io/)
* [QuPath Read the docs](https://qupath.readthedocs.io/en/stable/index.html)


Important Notes
This application generally requires a [graphical connection using NX](https://hpc.nih.gov/docs/nx.html)


* Module Name: QuPath (see [the modules page](/apps/modules.html) for more information)
 * Environment variable set 
	+ QUPATH\_ROOT* Example files in $QUPATH\_ROOT/TESTDATA



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Start a session in [NX](https://hpc.nih.gov/docs/nx.html), open a terminal, allocate an [interactive session](/docs/userguide.html#int), and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf ~]$ **sinteractive -c8 --mem=16g --gres=lscratch:10**
salloc: Pending job allocation 30957652
salloc: job 30957652 queued and waiting for resources
salloc: job 30957652 has been allocated resources
salloc: Granted job allocation 30957652
salloc: Waiting for resource configuration
salloc: Nodes cn0849 are ready for job

[user@cn0849 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn0849 30957652]$ **module load QuPath**
[+] Loading QuPath  0.3.2  on cn0849 
[+] Loading singularity  3.8.5-1  on cn0849 

[user@cn0849 30957652]$ **cp -r $QUPATH\_ROOT/TESTDATA .**

[user@cn0849 30957652]$ **cd TESTDATA/**

[user@cn0849 TESTDATA]$ **ls**
CMU-1-JP2K-33005.svs  CMU-1-Small-Region.svs  CMU-1.svs  CMU-2.svs  CMU-3.svs  index.yaml  JP2K-33003-1.svs  JP2K-33003-2.svs

[user@cn0849 TESTDATA]$ **QuPath.sh --image CMU-1-JP2K-33005.svs**
OpenJDK 64-Bit Server VM warning: Option --illegal-access is deprecated and will be removed in a future release.
Jan 25, 2022 9:04:02 PM com.sun.javafx.application.PlatformImpl startup
WARNING: Unsupported JavaFX configuration: classes were loaded from 'unnamed module @616fe72b'
21:04:04.013 [JavaFX Application Thread] [INFO ] qupath.lib.common.ThreadTools - Setting parallelism to 7
21:04:04.336 [JavaFX Application Thread] [INFO ] qupath.lib.gui.QuPathGUI - QuPath build: Version: 0.3.2
Build time: 2022-01-17, 08:49
Latest commit tag: '71884c6'
21:04:04.337 [JavaFX Application Thread] [INFO ] qupath.lib.gui.QuPathGUI - Setting tile cache size to 48312.00 MB (25.0% max memory)

(QuPath:11973): Gdk-WARNING **: 16:04:05.273: XSetErrorHandler() called with a GDK error trap pushed. Don't do that.
21:04:05.888 [JavaFX Application Thread] [INFO ] qupath.lib.scripting.QP - Initializing type adapters
21:04:06.648 [JavaFX Application Thread] [INFO ] q.l.i.s.b.BioFormatsOptionsExtension - Bio-Formats version 6.7.0
21:04:06.651 [JavaFX Application Thread] [INFO ] qupath.lib.gui.QuPathGUI - Loaded extension Bio-Formats options (Bio-Formats 6.7.0) (18 ms)
21:04:06.721 [JavaFX Application Thread] [INFO ] qupath.lib.gui.QuPathGUI - Loaded extension ImageJ extension (69 ms)
21:04:06.750 [JavaFX Application Thread] [INFO ] qupath.lib.gui.QuPathGUI - Loaded extension Processing extension (29 ms)
21:04:06.830 [JavaFX Application Thread] [INFO ] qupath.lib.gui.QuPathGUI - Loaded extension Rich script editor extension (79 ms)
21:04:06.831 [JavaFX Application Thread] [INFO ] qupath.lib.gui.QuPathGUI - Loaded extension SVG export extension (1 ms)
21:04:06.877 [JavaFX Application Thread] [INFO ] q.l.i.s.o.OpenslideServerBuilder - OpenSlide version 3.4.1
21:04:07.549 [JavaFX Application Thread] [INFO ] qupath.lib.gui.QuPathApp - Starting QuPath with parameters: [--image=CMU-1-JP2K-33005.svs]
21:04:07.550 [qupathgui-1] [INFO ] qupath.lib.gui.QuPathGUI - Update check for https://github.com/qupath/qupath
21:04:07.910 [qupathgui-1] [ERROR] qupath.lib.gui.QuPathGUI - Update check failed for QuPath (owner=qupath, repo=qupath)
21:04:07.965 [JavaFX Application Thread] [INFO ] q.l.i.s.b.BioFormatsServerOptions - Setting max Bio-Formats readers to 8
21:04:08.245 [JavaFX Application Thread] [WARN ] q.l.i.s.b.BioFormatsImageServer$ReaderPool - Temp memoization directory created at /tmp/qupath-memo-4974945218509452371
21:04:08.245 [JavaFX Application Thread] [WARN ] q.l.i.s.b.BioFormatsImageServer$ReaderPool - If you want to avoid this warning, either specify a memoization directory in the preferences or turn off memoization by setting the time to < 0
21:04:12.187 [JavaFX Application Thread] [INFO ] qupath.lib.gui.viewer.QuPathViewer - Image data set to ImageData: Not set, CMU-1-JP2K-33005.svs

```

  

At this point you should see an Graphical User Interface like the one pictured below:
  

  


![QuPath GUI image](/images/QuPath.jpg)


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
In addition to running through the GUI, QuPath.sh can be used to run scripts non-interactively. Create a groovy script with the image processing commands you want to run e.g. myScript. Then create a batch input file (e.g. QuPath\_batch.sh). For example:



```

#!/bin/bash
set -e
module load QuPath
QuPath.sh script myScript

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] QuPath_batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. QuPath.swarm). For example:



```

QuPath.sh script Script1
QuPath.sh script Script2
QuPath.sh script Script3
QuPath.sh script Script4

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f QuPath.swarm [-g #] [-t #] --module QuPath
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module QuPath Loads the QuPath module for each subjob in the swarm 
 | |
 | |
 | |








