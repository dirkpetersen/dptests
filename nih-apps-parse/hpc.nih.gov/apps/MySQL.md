

document.querySelector('title').textContent = 'MySQL';
MySQL


|  |
| --- |
| 
Quick Links
[Mysql database on biobase](#biobase)
[Temporary User-Managed Mysql](#temp)
[Permanent User-Managed Mysql](#perm)
 |


[MySQL](http://www.mysql.com/) is a general-purpose RDBMS, and is widely used as a backend for many applications. On Biowulf, users can either request a database on the mysql host biobase.nih.gov, or run a temporary Mysql instance on a node, or run an occasional Mysql instance pointing to a permanent location in the user's data area. These options are detailed below.




Choose your modality


|  |  |  |
| --- | --- | --- |
| **Modality** | **Pros** | **Cons** |
| [biobase.nih.gov](#biobase) | * always available
 * accessible within NIH network
 * most reliable data storage
 | * no root access
 * strict limits on databases and use
 * controlled and restricted by HPC Staff
 |
| [local\_mysql](#temp) | * easy to create and manage
 * full root access
 * full control over administration
 * fastest and most reliable access
 | * must be started within an HPC job
 * not accessible outside cluster
 * temporary instance on /lscratch
 * data must be archived between jobs
 |
| [local\_mysql --basedir](#perm) | * easy to create and manage
 * full root access
 * full control over administration
 * no need to archive data
 | * must be started within an HPC job
 * not accessible outside cluster
 * data stability not 100% reliable
* access speed affected by network slowdowns
* InnoDB engine not available
 |



Mysql database on biobase.nih.gov
The mysql host biobase.nih.gov is available for user accounts. This server is available only
within the NIH network, and is intended only for use with applications run on HPC @ NIH hosts.
To obtain a mysql account on biobase.nih.gov, please contact [staff@hpc.nih.gov](mailto:staff@hpc.nih.gov).


One of the primary purposes for biobase.nih.gov is hosting the backend databases for our local
mirror of the [UCSC Genome Browser](http://genome.cit.nih.gov/). Users can access the
local databases using a special read-only account. Please contact [staff@hpc.nih.gov](mailto:staff@hpc.nih.gov) for details.


#### Using biobase from the cluster


Because nodes in the cluster are sequestered to an internal network, the host 'biobase.nih.gov' should
not be used as the mysql host. Instead, use 'biobase':



```
[node1]$ mysql -u *user* -p -h **biobase**
```





Temporary Mysql instance
Once a mysql module has been loaded, the command **local\_mysql** allows a mysql
instance to be started on a node. The **local\_mysql** has four sub-commands:


* **create** -- creates the MySQL installation tree
 * **start** -- starts the server
 * **stop** --stops the server
 * **destroy** -- removes the MySQL installation tree

 * **archive** -- create a tarball of basedir
 * **restore** -- restore basedir from tarball
 * **status** -- report status of local mysql server


Here is an example of how to start a mysql server on an interactive
node, using the dynamically-allocated local /lscratch disk as the root directory.


#### Start interactive session


In this example, we allocate 10 GB of local scratch space for the session. 
See <https://hpc.nih.gov/docs/userguide.html#local> for more information about using /lscratch.
See <https://hpc.nih.gov/docs/userguide.html#int> for more information on allocating interactive sessions.
*Job and node ids have been changed to protect the innocent*.




```

[biowulf]$ **sinteractive --mem=10g --gres=lscratch:10**
salloc.exe: Pending job allocation ***123456789***
salloc.exe: job ***123456789*** queued and waiting for resources
salloc.exe: job ***123456789*** has been allocated resources
salloc.exe: Granted job allocation ***123456789***
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes ***node1*** are ready for job

```

#### Initialize and start a mysql instance


Once the interactive session is ready, load the mysql module of choice, initialize and start the mysql instance.


By default, **local\_mysql** will set the root password of the instance to
mysql123. This can be changed using the option **--password**.



```

[node1]$ **module load mysql**
[node1]$ **local\_mysql create**
```

The default engine is **MyISAM**. This can be changed by editing **my.cnf** prior to starting the server:



```

[node1]$ **cat /lscratch/$SLURM\_JOB\_ID/mysql/my.cnf**
...
default_storage_engine      = *MyISAM*
default_tmp_storage_engine  = *MyISAM*
...
[node1]$ **vi /lscratch/$SLURM\_JOB\_ID/mysql/my.cnf**
[node1]$ **cat /lscratch/$SLURM\_JOB\_ID/mysql/my.cnf**
...
default_storage_engine      = *InnoDB*
default_tmp_storage_engine  = *InnoDB*
...

```

Set **MYSQL\_HOME** and **MYSQL\_UNIX\_PORT** variables:



```
[node1]$ **export MYSQL\_HOME=/lscratch/$SLURM\_JOB\_ID/mysql**
[node1]$ **export MYSQL\_UNIX\_PORT=/lscratch/$SLURM\_JOB\_ID/mysql/mysql.sock**
```

Now start up the server:



```

[node1]$ **local\_mysql start**
MySQL is running now
You can now log in using 'mysql -u root -p'mysql123'

```

By default, the initialization procedure creates a new directory under /lscratch (**/lscratch/$SLURM\_JOB\_ID/mysql**) and fills it with files and directories:



```

[node1]$ **ls -F /lscratch/*123456789*/mysql**
bin/  data/  my.cnf  mysql.err  mysql.pid  mysql.sock  share/  support-files/  tmp/

```

* **bin/** - holds executables needed for mysql server
* **data/** - location of databases and tables
* **my.cnf** - personal configuration
* **mysql.err** - error logfile
* **mysql.pid** - holds process id for mysql server instance
* **mysql.sock** - socket file for mysql server instance
* **share/** - miscellaneous files
* **support-file/** - miscellaneous files
* **tmp/** - temporary files needed for mysql server


#### Log in do a few things


Once the mysql instance is up and running, you can log in using either:


* **The full socket file path:** --socket=/lscratch/$SLURM\_JOB\_ID/mysql/mysql.sock
* **The host and port:** --host=$SLURM\_NODELIST --port=55555


By default, **local\_mysql** will set the port to 55555. **This is different than the standard
mysql port of 3306 and therefore must be explicitly specified either by using the socket file or with the --port option.** 
The port can be changed using the option **--port**.



```

[node1]$ **mysql -u root -p'mysql123'**
Welcome to the MySQL monitor.  Commands end with ; or \g.
Your MySQL connection id is 3
Server version: 5.7.30 Source distribution

Copyright (c) 2000, 2020, Oracle and/or its affiliates. All rights reserved.

Oracle is a registered trademark of Oracle Corporation and/or its
affiliates. Other names may be trademarks of their respective
owners.

Type 'help;' or '\h' for help. Type '\c' to clear the current input statement.

mysql> **CREATE USER me IDENTIFIED BY 'pass123';**
mysql> **CREATE DATABASE `my\_db`;**
mysql> **GRANT ALL ON `my\_db`.\* TO me;**
mysql> **USE my\_db;**
mysql> **CREATE TABLE my\_table (x int, y int, z int);**
mysql> **DESC my\_table;**
+-------+---------+------+-----+---------+-------+
| Field | Type    | Null | Key | Default | Extra |
+-------+---------+------+-----+---------+-------+
| x     | int(11) | YES  |     | NULL    |       |
| y     | int(11) | YES  |     | NULL    |       |
| z     | int(11) | YES  |     | NULL    |       |
+-------+---------+------+-----+---------+-------+
mysql> **exit**
[node1]$

```

#### Simplify using ~/.my.cnf


MySQL can read default information from configuration files. For personal use, the default is **~/.my.cnf**. This file can be edited
to simplify use of mysql.



```

# Configuration for local_mysql
[mysql]
auto-rehash=FALSE
prompt= \h:\d>\_

[client]
user=me
password=pass123
database=my_db
host=node1
port=55555
socket=/lscratch/123456789/mysql/mysql.sock

```

Note that the **host**, **port** and **socket** values will likely be different.


Creating this personal **~/.my.cnf** file will allow simple connections:



```
[node1]$ mysql
...
localhost:my_db>

```

#### Log in from another machine


While the interactive session is running, the mysql instance will be accessible from any node within the Biowulf cluster
using the **--host** and **--port** options.



```
[node2]$ **mysql -u me -p'pass123' --host=*node1* --port=55555**
```

By default the root account initially created **does not allow access from outside the original host where mysql was initialized**.
Either modify the root account (not good idea) or create an alternate account (better idea) for access from other nodes.


#### Stop mysql prior to archiving and exiting the interactive session



```

[node1]$ **local\_mysql stop**

MySQL is stopped!

```

#### Archive the contents


The contents of the temporary server instance can be archived into a single tarball for restarting later and on a different node. This is the recommended way of keeping instances long-term.



```

[node1]$ **local\_mysql archive --archivefile=/data/$USER/save.tgz**
[node1]$ **ls -l /data/$USER/save.tgz**
-rw-r--r-- 1 user user 1371882 Jul  1 11:35 /data/user/save.tgz

```

#### Exit the interactive session


Lastly, exit the interactive session. The contents of /lscratch will be destroyed!



```

[node1]$ **exit**
salloc.exe: Relinquishing job allocation *123456789*
[biowulf]$

```

#### Start a new interactive session, restore the contents, and start mysql again



```

[biowulf]$ **sinteractive --mem=10g --gres=lscratch:10**
salloc.exe: Pending job allocation ***234567891***
... 
[node3]$ **module load mysql**
[node3]$ **local\_mysql restore --archivefile=/data/$USER/save.tgz**
[node3]$ **local\_mysql start**

MySQL is running now
You can now log in using 'mysql -u root -p'mysql123' --socket=/lscratch/***234567891***/mysql/mysql.sock'

```




Mysql database with permanent location in /data
If the mysql files are required permanently, such that the server could be brought up again
and again, **and archiving is not a possibility (see above)**, then don't take the default **--basedir** value. Instead, point the **--basedir** option to a directory in 
a /data area.


By default mysql files are written to /lscratch. This can be changed using the option **--basedir**.
**local\_mysql** will attempt to create the directory path given by **--basedir**, and will complain if
the path already exists.



```

[node1]$ **module load mysql**
[node1]$ **local\_mysql --basedir /data/*user*/mysql start**
...
[node1]$ **mysql -u root -p'mysql123' --host=$SLURM\_NODELIST --port=55555**
mysql>

```

And as always be sure to stop the server when you're done using it:



```
[node1]$ **local\_mysql --basedir /data/*user*/mysql stop**
```

Here is a script that can be used to launch a mysql server in a batch job. It traps signals so that the server is cleanly shut down prior to the job ending.



```

#!/bin/bash
#SBATCH --time=8:00:00 --signal=B:USR1@60 --job-name=mysql_server

term_handler()
{
    [[ -f /data/*user*/mysql/mysql.pid ]] && local_mysql --basedir /data/*user*/mysql --port 55555 stop
    exit
}

# associate the function "term_handler" with signals
trap 'term_handler' TERM
trap 'term_handler' EXIT
trap 'term_handler' USR1

# start the server
module load mysql
local_mysql --basedir /data/*user*/mysql --port 55555 start

# sleep for the time with sbatch (above)
sleep 8h & wait
[[ -f /data/*user*/mysql/mysql.pid ]] && local_mysql --basedir /data/*user*/mysql --port 55555 stop

```

To find the node on which the mysql server is running, use the jobname as a filter:



```
dashboard_cli jobs --jobname mysql_server --fields nodelist --raw --no-header
```


#### Restart the mysql server


Then later on, the server can be started up with the original database files:



```
[biowulf]$ **sinteractive**
...
[node1]$ **module load mysql**
[node1]$ **local\_mysql --basedir /data/*user*/mysql start**

```

For more information about how to use **local\_mysql**, type



```
local_mysql --help
```







