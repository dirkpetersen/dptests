

document.querySelector('title').textContent = 'Apache Spark on Biowulf';
Apache Spark on Biowulf
Apache Spark is a large-scale data processing engine that performs
in-memory computing. Spark offers bindings in Java, Scala, Python and R for
building parallel applications.



Documentation
* Apache Spark [home page](https://spark.apache.org/)
* [Documentation](https://spark.apache.org/docs/latest/)



Manage Spark clusters on biowulf
Spark clusters running on biowulf nodes can be managed with the `spark`
tool. Once a Spark cluster has been started, it can be used interactively or it can
be used to submit Spark jobs to. Once there is no more need for the cluster
it must be shut down.



```

[user@biowulf]$ **module load spark**
[user@biowulf]$ **spark**

NAME
        spark - administer spark clusters on compute nodes

SYNOPSIS
        spark cmd [options]

COMMANDS
        help   - show this help or help for specific cmd
        start  - start a new spark cluster
        list   - list clusters
        show   - show cluster details
        stop   - shut down a cluster
        clean  - clean up the directory where cluster
                 info is stored (/data/user/.spark-on-biowulf)

DESCRIPTION
        This tool is used to start, stop, monitor, and
        use spark clusters running on compute nodes.

[user@biowulf]$ **spark help start**

NAME
      spark start - start new spark cluster

SYNOPSIS
      spark start [options] nnodes

DESCRIPTION
      Provisions a new spark cluster with  nodes

 -t M max runtime for the cluster in minutes. Minimum is
 10. Default is 30. Actual runtime of the cluster
 is slightly less to allow for startup and clean
 shutdown
 -l copy the spark node logs back to the spark cluster
 directory in shared space

```

Let's start a spark cluster on 2 nodes. This tool uses 56 CPU nodes
with 256GB of memory and set a max runtime of 2h.



```

[user@biowulf]$ **spark start -t 120 2**
INFO: Submitted job for cluster TkJvrN

```

The `spark` tool stores information about all its clusters
in `/data/$USER/.spark-on-biowulf`. That includes clusters
that already completed.



```

[user@biowulf]$ **tree /data/$USER/.spark-on-biowulf**
/data/user/.spark-on-biowulf
`-- [user   4.0K]  TkJvrN
    |-- [user   4.1K]  jobscript.sh
    |-- [user   4.0K]  logs
    |-- [user    109]  prop
    `-- [user      9]  slurm_job_id

```

We can check on the status of our clusters with



```

[user@biowulf]$ **spark list**
Cluster id  Slurm jobid                state
---------- ------------ --------------------
    TkJvrN     18256246              RUNNING

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


```

Note that it may take a couple of minutes after the cluster starts running for
the detailed information to be complete. Now that the cluster is running, we can
put it to work. For example, lets use the cluster interactively with pyspark and
have a bit of a look at the sqlite source code:



```

[user@biowulf]$ module load python/3.8
[user@biowulf]$ pyspark --master spark://cn0443:7077 --executor-memory=40g
Python 3.8.12 | packaged by conda-forge | (main, Mar 24 2022, 23:25:59)
[GCC 10.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
[...some warnings...]
Welcome to
      ____              __
     / __/__  ___ _____/ /__
    _\ \/ _ \/ _ `/ __/  '_/
   /__ / .__/\_,_/_/ /_/\_\   version 3.2.1
      /_/

Using Python version 3.8.12 (main, Mar 24 2022 23:25:59)
Spark context Web UI available at http://biowulf2-e0:4040
Spark context available as 'sc' (master = spark://cn4170:43701, app id = app-20220518160041-0000).
SparkSession available as 'spark'.

>>> txt = spark.sparkContext.textFile("sqlite3.c")
>>> txt.count()
202918
>>> defines = txt.filter(lambda l: l.startswith('#define'))
>>> defines.count()
2492
>>> defines.first()
u'#define SQLITE_CORE 1'
>>> txt.map(lambda l: len(l)).reduce(lambda a, b: a if (a>b) else b)
260
>>> Ctrl-D

[user@biowulf]$

```

pyspark will use the first python interpreter on the path. That means
loading the python/3.8 module, for example, will give you a python 3.8 spark
shell. There are similar shells for R (`sparkR`) and scala
(`spark-shell`).


pyspark can also be used with a jupyter notebook:



```

[user@biowulf]$ **sinteractive --mem=10g --tunnel**
...
[user@cn3144]$ **module load spark jupyter**
[user@cn3144]$ **pyspark-jupyter --master spark://cn0443:7077 --executor-memory=40g**
Running on port 10762 on cn3421 listening on localhost
[...snip...]

```

See our [Jupyter documentation](https://hpc.nih.gov/apps/jupyter.html)
for information about how to connect to a jupyter notebook running on a compute
node.


And we can uses `spark-submit` to submit spark jobs to the
cluster



```

 [user@biowulf]$ **spark-submit \
 --driver-memory=3g \
 --master spark://cn0443:7077 \
 --deploy-mode client \
 --executor-cores=2 \
 --executor-memory=3g \
 ./pi.py**

[...snip...]
/pi.py:36, took 562.603120 s
PI=3.141584232000

```

`spark clean` will delete metadata of finished spark clusters
from the spark metadata directory.



Local pseudo clusters for debugging
For development and debugging it can be convenient to run spark in local
pseudo cluster mode. This is how the interactive shells start up when
no master is provided on the command line. In order to use this a node
has to be allocated exclusively.



```

[user@biowulf ~]$ **sinteractive --exclusive --ntasks=1 --cpus-per-task=32 \
 --constraint cpu32**
[user@cn0182 ~]$ **module load spark**
[user@cn0182 ~]$ **pyspark**

```

After some initialization output, you will see the following:



```

Welcome to
      ____              __
     / __/__  ___ _____/ /__
    _\ \/ _ \/ _ `/ __/  '_/
   /__ / .__/\_,_/_/ /_/\_\   version 2.4.0
      /_/

Using Python version 3.6.8 (default, Dec 30 2018 01:22:34)
SparkSession available as 'spark'.
>>> spark.sparkContext.master
'local[32]'
>>>

```



