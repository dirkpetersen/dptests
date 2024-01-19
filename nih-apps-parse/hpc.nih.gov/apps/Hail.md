

document.querySelector('title').textContent = 'Hail on Biowulf';
Hail on Biowulf


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



Hail is an open-source, scalable framework for exploring and analyzing genomic data. See <https://hail.is/docs/0.2/index.html> for more information. 



Documentation
* [Hail Main Site](https://hail.is/docs/0.2/tutorials-landing.html)


Important Notes
* Module Name: hail (see [the modules page](/apps/modules.html) for more information)
* Cluster/distibuted computing



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive -c 16 --mem 40g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ]$ **module load hail**
[+] Loading hail  0.2.3 on cn3344 
[+] Loading singularity  on cn3344 

[user@cn3144]$ **wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/\
20100804/ALL.2of4intersection.20100804.sites.vcf.gz**
[user@cn3144 ]$ **ipython**
Python 3.6.7 (default, Oct 25 2018, 09:16:13) 
Type 'copyright', 'credits' or 'license' for more information
IPython 7.1.1 -- An enhanced Interactive Python. Type '?' for help.

In [1]: **import hail as hl**

In [2]: **hl.init()**                                                                                            
using hail jar at /usr/local/lib/python3.6/dist-packages/hail/hail-all-spark.jar
Using Spark's default log4j profile: org/apache/spark/log4j-defaults.properties
Setting default log level to "WARN".
To adjust logging level use sc.setLogLevel(newLevel). For SparkR, use setLogLevel(newLevel).
Running on Apache Spark version 2.2.2
SparkUI available at http://10.2.9.172:4040
Welcome to
     __  __     <>__
    / /_/ /__  __/ /
   / __  / _ `/ / /
  /_/ /_/\_,_/_/_/   version 0.2-a2eaf89baa0c
LOGGING: writing to /spin1/scratch/teacher/hail-20181129-1008-0.2-a2eaf89baa0c.log

In [4]: **hl.import\_vcf('ALL.2of4intersection.20100804.sites.vcf.gz',force\_bgz=True).write('sample.vds')**                                                                                  
[Stage 1:===============================>                          (7 + 6) / 13]2018-11-29 15:10:40 Hail: INFO: Coerced sorted dataset
[Stage 2:===========================================>             (10 + 3) / 13]2018-11-29 15:11:47 Hail: INFO: wrote 25488488 items in 13 partitions to sample.vds

[user@cn3144 ]$ **exit**salloc.exe: Relinquishing job allocation 46116226

```

 Run hail with jupyter notebook on single node:

```

[user@biowulf ]$**sinteractive -c 16 --mem 40g --tunnel**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job
Created 1 generic SSH tunnel(s) from this compute node to                  
biowulf for your use at port numbers defined                               
in the $PORTn ($PORT1, ...) environment variables.                         
                                                                           
                                                                           
Please create a SSH tunnel from your workstation to these ports on biowulf.
On Linux/MacOS, open a terminal and run:                                   
                                                                           
    ssh  -L 33327:localhost:33327 biowulf.nih.gov                          
                                                                           
For Windows instructions, see https://hpc.nih.gov/docs/tunneling          
[user@cn3144]$ **module load hail**
[user@cn3144]$ **wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/\
20100804/ALL.2of4intersection.20100804.sites.vcf.gz**
[user@cn3144]$ **jupyter notebook --ip localhost --port $PORT1 --no-browser**
[I 17:11:40.505 NotebookApp] Serving notebooks from local directory
[I 17:11:40.505 NotebookApp] Jupyter Notebook 6.4.10 is running at:
[I 17:11:40.505 NotebookApp] http://localhost:37859/?token=xxxxxxxx
[I 17:11:40.506 NotebookApp]  or http://127.0.0.1:37859/?token=xxxxxxx
[I 17:11:40.506 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[C 17:11:40.512 NotebookApp]

    To access the notebook, open this file in a browser:
        file:///home/apptest1/.local/share/jupyter/runtime/nbserver-29841-open.html
    Or copy and paste one of these URLs:
        http://localhost:37859/?token=xxxxxxx
     or http://127.0.0.1:37859/?token=xxxxxxx


```

 Then you can open a browser from your computer to connect to the jupyter notebook:
 ![hail_jupyter](/images/hail_jupyter.png)


 Run hail with jupyter notebook with spark(3.1.3):

```

[user@biowulf]$ **module load spark/3.1.3**
[user@biowulf]$ **spark start -t 120 2**
INFO: Submitted job for cluster TkJvrN
[user@biowulf]$ **spark list -d**
Cluster id  Slurm jobid                state
---------- ------------ --------------------
    TkJvrN     18256246              RUNNING
               nodes: 2
            max_time: 120
               spark: 2.4.0
              job_id: 18256246
               start: 2019-01-16 12:33:03
            nodelist: cn[3769-3770]
              master: spark://cn3769:7077
        master_webui: http://cn3769:8080
              tunnel: ssh -L 8080:cn3769:8080 -N
[user@biowulf]$**sinteractive --tunnel**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job
Created 1 generic SSH tunnel(s) from this compute node to
biowulf for your use at port numbers defined
in the $PORTn ($PORT1, ...) environment variables.


Please create a SSH tunnel from your workstation to these ports on biowulf.
On Linux/MacOS, open a terminal and run:

    ssh  -L 33327:localhost:33327 biowulf.nih.gov

For Windows instructions, see https://hpc.nih.gov/docs/tunneling
[user@cn3144]$ **module load hail/0.2.95**
[user@cn3144]$ **hail-jupyter --master spark://cn3769:7077 \
 --driver-memory=6g \
 --driver-cores=2 \
 --executor-cores=2 \
 --num-executors=50 \
 --executor-memory=5g**
[I 17:11:40.505 NotebookApp] Serving notebooks from local directory
[I 17:11:40.505 NotebookApp] Jupyter Notebook 6.4.10 is running at:
[I 17:11:40.505 NotebookApp] http://localhost:37859/?token=xxxxxxxx
[I 17:11:40.506 NotebookApp]  or http://127.0.0.1:37859/?token=xxxxxxx
[I 17:11:40.506 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[C 17:11:40.512 NotebookApp]

    To access the notebook, open this file in a browser:
        file:///home/apptest1/.local/share/jupyter/runtime/nbserver-29841-open.html
    Or copy and paste one of these URLs:
        http://localhost:33327/?token=xxxxxxx
     or http://127.0.0.1:33327/?token=xxxxxxx


```

 Then you can open a browser from your computer to connect to the jupyter notebook to run on a spark cluster:
 ![hail_jupyter](/images/hail_jupyter_spark.png)



Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a python script (e.g. hail-script.py). For example:



```

#!/usr/bin/env python3
import hail as hl
mt = hl.balding_nichols_model(n_populations=3,
                              n_samples=500,
                              n_variants=500_000,
                              n_partitions=32)
mt = mt.annotate_cols(drinks_coffee = hl.rand_bool(0.33))
gwas = hl.linear_regression_rows(y=mt.drinks_coffee,
                                 x=mt.GT.n_alt_alleles(),
                                 covariates=[1.0])
gwas.order_by(gwas.p_value).show(25)

```

Create a batch input file (e.g. hail.sh). For example:



```

#!/bin/bash
module load hail
python3-hail hail-script.py 

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch -c 16 --mem 40g hail.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. hail.swarm). For example:



```

hail1.py
hail2.py
hail3.py

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f hail.swarm -g 30 -t 16 --module hail
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module hail Loads the hail module for each subjob in the swarm 
 | |
 | |
 | |
















