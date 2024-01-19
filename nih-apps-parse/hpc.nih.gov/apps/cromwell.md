

document.querySelector('title').textContent = 'cromwell on Biowulf';
cromwell on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



A Workflow Management System geared towards scientific workflows.



Documentation
Cromwell is a complex tool for running workflows described by the workflow description
language (WDL). On Biowulf, cromwell is used to either run local workflows (all on
one node) or distributed workflows (each task a slurm job). Server mode is not
supported. Comprehensive documentation of cromwell or WDL is bejond the scope of
this brief document. See the following links for detailed documentation:


* [Cromwell](https://github.com/broadinstitute/cromwell) on GitHub
* [WDL home](https://software.broadinstitute.org/wdl/)
* [WDL Tutorials](https://software.broadinstitute.org/wdl/documentation/topic?name=wdl-tutorials)
* [WDL](https://github.com/broadinstitute/wdl) on GitHub


Important Notes
* Module Name: cromwell (see [the modules page](/apps/modules.html) for more information)
* Cromwell can run in local mode or it can use slurm to run tasks
* Important environment variables: `**$CROMWELL\_JAR**`,
 `**$CROMWELL\_TEST\_DATA**` , `**$WDLTOOL\_JAR**`, `**$WOMTOOL\_JAR**`
`**$CROMWELL\_CONFIG**`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive -c6 --mem=20g --gres=lscratch:10**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load cromwell**
[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **cp -r ${CROMWELL\_TEST\_DATA:-none}/\* .**
[user@cn3144]$ **tree -s**
.
|-- [user   4.0K]  data
|   |-- [user    15M]  sample1.fq.gz
|   |-- [user    15M]  sample2.fq.gz
|   |-- [user    15M]  sample3.fq.gz
|   |-- [user    15M]  sample4.fq.gz
|   |-- [user    15M]  sample5.fq.gz
|   `-- [user    15M]  sample6.fq.gz
|-- [user    189]  input.json
`-- [user   1.1K]  wf.wdl

```

First, let's look at the workflow



```

[user@cn3144]$ **cat wf.wdl**

task lengths {
  File fq_file
  String prefix = sub(fq_file, ".fq.gz", "")
  command <<<
    module load gnuplot
    zcat ${fq_file} \
        | awk 'NR % 4 == 2 {print length($1)}' \
        | sort \
        | uniq -c \
        | sort -k2,2n \
        | awk '{print $2"\t"$1}' \
        > ${prefix}.len
    gnuplot -e "set term dumb; set style data lines; plot '${prefix}.len'" > ${prefix}.len.plot
  >>>
  output {
    File dat  = "${prefix}.len"
    File plot = "${prefix}.len.plot"
  }
  runtime { rt_mem: 2000 rt_time: 10 }
}

task reads {
  File fq_file
  String prefix = sub(fq_file, ".fq.gz", "")
  command <<<
    zcat ${fq_file} \
        | awk 'NR % 4 == 2' \
        | sort \
        | uniq -c \
        | sort -k1,1nr \
        > ${prefix}.fq.nr
  >>>
  output {
    File dat  = "${prefix}.fq.nr"
  }
  runtime { rt_mem: 2000 rt_time: 10 }
}

workflow wf {
  Array[File] fq_files
  scatter(s in fq_files) {
    call lengths { input: fq_file=s }
    call reads {input: fq_file=s }
  }
}


```

Two tasks that summarize the length trimmed reads in input files, create an
ascii graph, and summarize all non-redundant reads. These are steps you might do
in a miRNA analysis. Let's validate the WDL file with womtool:



```

[user@cn3144]$ **java -jar ${WOMTOOL\_JAR} validate wf.wdl**
Success!
[user@cn3144]$ **echo $?**
0

```

So womtool returned a 0 and did not report any errors - the workflow is valid. Next
we can combined the workflow with a json file specifying the input files. To do that,
we need to modify the input.json file to reflext the current path:



```

[user@cn3144]$ **sed -i "s:XXX:$PWD:" input.json**
[user@cn3144]$ **cat input.json**
{
  "wf.fq_files": [
      "/path/to/here/data/sample1.fq.gz",
      "/path/to/here/data/sample2.fq.gz",
      "/path/to/here/data/sample3.fq.gz",
      "/path/to/here/data/sample4.fq.gz",
      "/path/to/here/data/sample5.fq.gz",
      "/path/to/here/data/sample6.fq.gz"
   ]

}
[user@cn3144]$ **java -jar ${CROMWELL\_JAR} run -i input.json wf.wdl**
[2017-09-08 16:04:32,62] [info] Slf4jLogger started
[2017-09-08 16:04:32,71] [info] RUN sub-command
[2017-09-08 16:04:32,71] [info]   WDL file: /path/to/here/wf.wdl
[2017-09-08 16:04:32,71] [info]   Inputs: /path/to/here/input.json
[...snip...]

[user@cn3144]$ **cat data/sample1.len.plot**
    400000 ++---------+----------+----------+---------+----------+---------++
           +'/path/to/here/data/sample1.len' ******                         +
    350000 ++             *                                                ++
           |              **                                                |
           |             * *                                                |
    300000 ++            * *                                               ++
           |             * *                                                |
    250000 ++            *  *                                              ++
           |             *  *                                               |
           |             *  *                                               |
    200000 ++           *   *                                              ++
           |            *    *                                              |
    150000 ++           *    *                                             ++
           |            *     *                                             |
           |           *      *                                             |
    100000 ++          *      *                                            ++
           |           *       *                                            |
     50000 ++         *        *                                           ++
           |         **         *                                           |
           +     **** +          **         +         +          +          +
         0 ******-----+----------+-*********************************-------++
           15         20         25         30        35         40         45


```

Similarly, the same workflow can be run by submitting tasks to slurm. For this, a
configuration file is needed. A working example can be found in ${CROMWELL\_CONFIG}.
It allows


* submission of regular and containerized jobs to Slurm
* Can specify `gpuCount` and `gpuType` in 
 the workflow and slurm will properly allocate the job and use the
 `--nv` option for singularity.
* configures the file system to use hard links, cached copies, or
 copies for localization of files
* uses call caching with a locally persisted Hsqldb database


For this simple workflow that file will suffice. For a more complicated workflow
you may have to copy it and modify according to your needs. Some potential changes
may include using a MySQL database backend or modifying the runtime property names
used in the configuration and raising the concurrent job limit from a low 10 jobs
in this simple example.


Note that in this mode, the current
 working directory has to be on a shared file system since information is
exchanged between jobs and the main cromwell process.


The singularity module has to be loaded before running a job that uses docker or
singularity containers.



```

[user@cn3144]$ **cd /data/$USER**
[user@cn3144]$ **mkdir cromwell\_test**
[user@cn3144]$ **cd cromwell\_test**
[user@cn3144]$ **cp -r ${CROMWELL\_TEST\_DATA:-none}/\* .**
[user@cn3144]$ **cat ${CROMWELL\_CONFIG}**
# include the application.conf at the top
include required(classpath("application"))

system {
  job-rate-control {
    jobs = 1
    per = 1 second
  }
}

database {
  profile = "slick.jdbc.HsqldbProfile$"
  db {
    driver = "org.hsqldb.jdbcDriver"
    url = """
    jdbc:hsqldb:file:cromwell-executions/cromwell-db/cromwell-db;
    shutdown=false;
    hsqldb.default_table_type=cached;hsqldb.tx=mvcc;
    hsqldb.result_max_memory_rows=10000;
    hsqldb.large_data=true;
    hsqldb.applog=0;
    hsqldb.lob_compressed=true;
    hsqldb.script_format=3
    """
    connectionTimeout = 120000
    numThreads = 2
   }
}

call-caching {
  enabled = true
  invalidate-bad-cache-results = true
}

backend {
  default = "Slurm"
  providers {
    Slurm {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        concurrent-job-limit = 10
        # If an 'exit-code-timeout-seconds' value is specified:
        #     - check-alive will be run at this interval for every job
        #     - if a job is found to be not alive, and no RC file appears after this interval
        #     - Then it will be marked as Failed.
        ## Warning: If set, Cromwell will run 'check-alive' for every job at this interval
        exit-code-timeout-seconds = 360 
        filesystems {
         local {
           localization: [
             # soft link does not work for docker with --contain. Hard links won't work
             # across file systems
             "hard-link", "cached-copy", "copy"
           ]
         }
        }
        default-runtime-attributes {
            maxRetries = 0
        }

        runtime-attributes = """
        Int runtime_minutes = 600
        Int cpu = 2
        # the _mb is meaningful and and can result in implicit conversions.
        Int memory_mb = 4000
        String queue = "norm"
        Int? gpuCount
        String? gpuType
        String? docker
        """

        submit = """
            sbatch \
              -J ${job_name} \
              -D ${cwd} \
              -o ${out} \
              -e ${err} \
              -t ${runtime_minutes} \
              -c ${cpu} \
              --mem ${memory_mb} \
              --partition ${queue} \
              ${if defined(gpuCount) then 
                        (if defined(gpuType) then ('--gres=gpu:' + gpuType + ':' + gpuCount)
                                             else ('--gres=gpu:' + gpuCount))
                        else ''} \
              --wrap "/bin/bash ${script}"
        """

        submit-docker = """
            # SINGULARITY_CACHEDIR needs to point to a directory accessible by
            # the jobs (i.e. not lscratch). Might want to use a workflow local
            # cache dir like in run.sh
            if [ -z $SINGULARITY_CACHEDIR ]; then
                CACHE_DIR=$HOME/.singularity
            else
                CACHE_DIR=$SINGULARITY_CACHEDIR
            fi
            mkdir -p $CACHE_DIR  
            LOCK_FILE=$CACHE_DIR/singularity_pull_flock

            # we want to avoid all the cromwell tasks hammering each other trying
            # to pull the container into the cache for the first time. flock works
            # on GPFS, netapp, and vast (of course only for processes on the same
            # machine which is the case here since we're pulling it in the master
            # process before submitting).
            flock --exclusive --timeout 1200 $LOCK_FILE \
                singularity exec --containall docker://${docker} \
                echo "successfully pulled ${docker}!" &> /dev/null

            sbatch \
              -J ${job_name} \
              -D ${cwd} \
              -o ${out} \
              -e ${err} \
              -t ${runtime_minutes} \
              -c ${cpu} \
              --mem ${memory_mb} \
              --partition ${queue} \
              ${if defined(gpuCount) then 
                        (if defined(gpuType) then ('--gres=gpu:' + gpuType + ':' + gpuCount)
                                             else ('--gres=gpu:' + gpuCount))
                        else ''} \
              --wrap "singularity exec ${if defined(gpuCount) then '--nv ' else ''} --containall --bind ${cwd}:${docker_cwd} docker://${docker} ${job_shell} ${docker_script}"
        """
        kill = "scancel ${job_id}"
        check-alive = "dashboard_cli jobs --is-active -j ${job_id} &> /dev/null"
        job-id-regex = "(\\d+)"
      }
    }
  }
}
[user@cn3144]$ **java -Dconfig.file=${CROMWELL\_CONFIG} \
 -jar ${CROMWELL\_JAR} run -i input.json wf.wdl**
[2017-09-08 16:08:25,56] [info] Slf4jLogger started
[2017-09-08 16:08:25,64] [info] RUN sub-command
[2017-09-08 16:08:25,64] [info]   WDL file: /spin1/users/user/test_data/cromwell/fnord/wf.wdl
[2017-09-08 16:08:25,64] [info]   Inputs: /spin1/users/user/test_data/cromwell/fnord/input.json
[...snip...]

```

End the sinteractive session



```

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. cromwell.sh). For example:



```

#!/bin/bash
module load cromwell
java -jar ${CROMWELL_JAR} run -i input.json wf.wdl

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] cromwell.sh
```







