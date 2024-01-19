

document.querySelector('title').textContent = 'nextflow on Biowulf';
nextflow on Biowulf



|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



Nextflow is a domain specific language modelled after UNIX pipes. It simplifies
writing parallel and scalable pipelines. The version installed on our systems
can run jobs locally (on the same machine) and by submitting to Slurm.


The code that is executed at each pipeline stage can be written in a number of
different languages (shell, python, R, ...).


Intermediate results for workflows are stored in the `$PWD/work`
directory which allows resuming execution of pipelines.


The language used to write pipeline scripts is an extension of 
[groovy](http://www.groovy-lang.org/).


Nextflow is a complex workflow management tool. Please read
the manual carefully and make sure to place appropriate limits on your pipeline
to avoid submitting too many jobs or running too many local processes.


Nextflow, when running many tasks appears to create many 
temp files in the ./work directory. Please make sure that your pipeline
does not inadvertently create millions of small files which would result
in a degradation of file system performance.


Common pitfalls
[top](#top)


**FATAL error while mount /gs6**
By 6.25.2023, /gs6 is retired from the cluster, thus if you see some errors like this: 
 
 WARNING: skipping mount of /gs6: stat /gs6: no such file or directory
 FATAL: container creation failed: mount /gs6->/gs6 error: while
 mounting /gs6: mount source /gs6 doesn't exist.
 
  

 Please update to the most recent version of nextflow.config, and then run your pipeline again:
 
 cp /usr/local/apps/nextflow/nextflow.config .
 
  


Documentation
* [Home page](http://www.nextflow.io/index.html)
* [Manual](http://www.nextflow.io/docs/latest/index.html)
* [GitHub](https://github.com/nextflow-io/nextflow)
* [slack](https://nextflow.slack.com)


Important Notes
* nf-core/2.9 (experimental) is packaged with nextflow/23.04.1, it is accessible after loading the newest nextflow module.

module load nextflow
 nf-core --help
 



* Module Name: nextflow (see [the modules page](/apps/modules.html) for more information)
* Nextflow can use local, slurm, or srun/ignite
* The master process submitting jobs should be run
 either as a batch job or on an interactive node - not on the biowulf
 login node.
* Please explicitly set the `pollInterval` and `queueStatInterval`
 to reduce the frequency
 with which nextflow polls slurm. The default frequency creates too many
 queries and results in unnecessary load on the scheduler.
* When your /home directory is full while running nextflow, it could be that the singularity cache is filling up, please redirect those to /data directory with:
 

 export NXF\_SINGULARITY\_CACHEDIR=/data/$USER/singularity;
 export SINGULARITY\_CACHEDIR=/data/$USER/.singularity;



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
First, let's do some basic local execution. For this we will allocate an interactive session:



```

[user@biowulf]$ **sinteractive --mem=10g -c2 --gres=lscratch:10**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job
[user@cn3144]$ **module load nextflow**

```

For the traditional hello world example we will parallelize the
uppercasing of different language greetings:



```

# create file of greetings
[user@cn3144]$ **mkdir testdir;cat > ./testdir/test1 <<EOF**
Hello world!
Hallo world!
Ciao world!
Salut world!
Bongiorno world!
Servus world!
**EOF**
[user@cn3144]$ **cat > ./testdir/test2 <<EOF**
Gruess Gott world!
Na was los world!
Gruetzi world!
Hello world!
Come va world!
Ca va world!
Hi world!
Good bye world!
**EOF**
[user@cn3144]$ **cat > test.R <<EOF**
args <- commandArgs(trailingOnly = TRUE)
library(readr)
df <- read_lines(args[1])

# create output
sink(args[2])

for (each in df){
  cat (toupper(each))
  cat ('\n')
  }

sink()
**EOF**

```

We then create a file called `rhello.nf` that describes
the workflow to be executed



```

// Declare syntax version
nextflow.enable.dsl=2

params.output_dir = './results'

process getsbatchlist {
  module 'R'

  input:
   path(input_file)
   each Rs

  publishDir "${params.output_dir}"

  output:
   path "${input_file}.txt"

  script:
   """
   Rscript ${Rs} ${input_file} ${input_file}.txt
   """
}

workflow {
  def inputf =  Channel.fromPath('./testdir/test*')
  def Rs = Channel.fromPath('./test.R')
  getsbatchlist(inputf,Rs) | view
}


```

The workflow is executed with



```

[user@cn3144]$ **nextflow rhello.nf**
N E X T F L O W  ~  version 23.04.1
Launching `test2.nf` [hopeful_cray] DSL2 - revision: 7401a333f4
executor >  local (2)
[28/6ccadb] process > getsbatchlist (2) [100%] 2 of 2 ✔
/gpfs/gsfs8/users/apptest2/work/82/9d153e5b2a5ab4399ab36beb01e552/test2.txt
/gpfs/gsfs8/users/apptest2/work/28/6ccadb35b01a2755c9670375ec1a05/test1.txt

[user@cn3144]$ **cat results/test1.txt**
HELLO WORLD!
HALLO WORLD!
CIAO WORLD!
SALUT WORLD!
BONGIORNO WORLD!
SERVUS WORLD!

```

Note that results are out of order.


The same workflow can be used to run each of the processes as a slurm job
by creating a `nextflow.config` file. We provide a file with correct
settings for biowulf at `/usr/local/apps/nextflow/nextflow.config`.
If you use this file please don't change settings for job submission and
querying (`pollInterval, queueStatInterval, and submitRateLimit`).
In particular you might want to remove the lscratch allocation if that does not apply to your workflow. Although
it was encouraged to use lscratch as much as you can.




```

[user@cn3144]$ **cp /usr/local/apps/nextflow/nextflow.config .**
[user@cn3144]$ **cat nextflow.config**

params {
  config_profile_description = 'Biowulf nf-core config'
  config_profile_contact = 'staff@hpc.nih.gov'
  config_profile_url = 'https://hpc.nih.gov/apps/nextflow.html'
  max_memory = '224 GB'
  max_cpus = 32
  max_time = '72 h'

}


// use a local executor for short jobs and it has to give -c and --mem to make nextflow
// allocate the resource automatically. For this the
// settings below may have to be adapted to the allocation for
// the main nextflow job.
executor {
    $local {
        queueSize = 100
        memory = "$SLURM_MEM_PER_NODE MB"
        cpus = "$SLURM_CPUS_PER_TASK"

    }
    $slurm {
        queue = 'norm'
        queueSize = 200
        pollInterval = '2 min'
        queueStatInterval = '5 min'
        submitRateLimit = '6/1min'
        retry.maxAttempts = 1
    }
}

singularity {
    enabled = true
    autoMounts = true
    cacheDir = "/data/$USER/singularity"
    envWhitelist='https_proxy,http_proxy,ftp_proxy,DISPLAY,SLURM_JOB_ID,SINGULARITY_BINDPATH'
}

env {
    SINGULARITY_CACHEDIR="/data/$USER/singularity"
    PYTHONNOUSERSITE = 1
}

profiles {
    biowulflocal {
        process.executor = 'local'
        process.cache = 'lenient'

    }

    biowulf {
        process {
            executor = 'slurm'
            maxRetries = 1

            clusterOptions = ' --gres=lscratch:200 '

            scratch = '/lscratch/$SLURM_JOB_ID'
            // with the default stageIn and stageOut settings using scratch can
            // result in humungous work folders
            // see https://github.com/nextflow-io/nextflow/issues/961 and
            //     https://www.nextflow.io/docs/latest/process.html?highlight=stageinmode
            stageInMode = 'symlink'
            stageOutMode = 'rsync'

            // for running pipeline on group sharing data directory, this can avoid inconsistent files timestamps
            cache = 'lenient'

        // example for setting different parameters for jobs with a 'gpu' label
        // withLabel:gpu {
        //    queue = 'gpu'
        //    time = '36h'
        //    clusterOptions = " --gres=lscratch:400,gpu:v100x:1 "
        //    containerOptions = " --nv "
        // }

        // example for setting different parameters for a process name
        //  withName: 'FASTP|MULTIQC' {
        //  cpus = 6
        //  queue = 'quick'
        //  memory = '6 GB'
        //  time = '4h'
        // }

        // example for setting different parameters for jobs with a resource label
        //  withLabel:process_low {
        //  cpus = 2
        //  memory = '12 GB'
        //  time = '4h'
        // }
        // withLabel:process_medium {
        //  cpus = 6
        //  memory = '36 GB'
        //  time = '12h'
        // }
        // withLabel:process_high {
        //  cpus = 12
        //  memory = '72 GB'
        //  time = '16 h'
        // }
     }
        timeline.enabled = true
        report.enabled = true
    }
}


[user@cn3144]$ **nextflow run -profile biowulf hello.nf**
N E X T F L O W  ~  version 20.10.0
Launching `hello.nf` [intergalactic_cray] - revision: f195027c60
executor >  slurm (15)
[34/d935ef] process > splitLetters        [100%] 1 of 1 ✔
HELLO WORLD!
[...snip...]
[97/85354f] process > convertToUpper (11) [100%] 14 of 14 ✔

```

Running nextflow with biowulf profile (slurm executor) using test input from nf-core:

```

[user@cn3144]$ **nextflow run nf-core/sarek -profile test,biowulf --outdir testout**
N E X T F L O W  ~  version 22.10.4
Launching `https://github.com/nf-core/sarek` [agitated_noyce] DSL2 - revision: c87f4eb694 [master]

WARN: Found unexpected parameters:
* --test_data_base: https://raw.githubusercontent.com/nf-core/test-datasets/modules
- Ignore this warning: params.schema_ignore_params = "test_data_base"



------------------------------------------------------
                                        ,--./,-.
        ___     __   __   __   ___     /,-._.--~'
  |\ | |__  __ /  ` /  \ |__) |__         }  {
  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                        `._,._,'
      ____
    .´ _  `.
   /  |\`-_ \      __        __   ___
  |   | \  `-|    |__`  /\  |__) |__  |__/
   \ |   \  /     .__| /¯¯\ |  \ |___ |  \
    `|____\´

  nf-core/sarek v3.1.2
------------------------------------------------------
...

```

Run nextflow with local executor (biowulflocal profile) to utilize allocated cpus and memory on compute node, --mem and -c is essential, and use lscratch as work directory for troubleshooting: 

```

[user@biowulf]$**sinteractive --mem=80g -c 32 --gres=lscratch:200**
[user@cn3144]$ **nextflow run nf-core/sarek
-r 3.2.3 \
-profile biowulflocal \
--wes \
--joint\_germline \
--input test.csv \
--tools haplotypecaller,vep,snpeff \
--outdir /data/$USER/sarek/ \
--genome GATK.GRCh38 \
--igenomes\_base /fdb/igenomes\_s3 \
--save\_output\_as\_bam \
-w /lscratch/$SLURM\_JOB\_ID \
--cache\_version 110 \
--vep\_cache /fdb/VEP/110/cache \
--snpeff\_cache /fdb/snpEff/5.1d/data/**

```

Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. sbatch\_nf\_main.sh) to run the master process. For example:



```

#! /bin/bash
#SBATCH --job-name=nextflow-main
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --gres=lscratch:200
#SBATCH --time=24:00:00

module load nextflow
export NXF_SINGULARITY_CACHEDIR=/data/$USER/singularity;
export SINGULARITY_CACHEDIR=/data/$USER/.singularity;
export TMPDIR=/lscratch/$SLURM_JOB_ID

nextflow run nf-core/sarek
-r 3.2.3 \
-profile biowulf \
--wes \
--joint_germline \
--input test.csv \
--tools haplotypecaller,vep,snpeff \
--outdir /data/$USER/sarek \
--genome GATK.GRCh38 \
--igenomes_base /fdb/igenomes_s3 \
--save_output_as_bam \
--cache_version 110 \
--vep_cache /fdb/VEP/110/cache/ \
--snpeff_cache /fdb/snpEff/5.1d/data/


```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch_nf_main.sh
```

The master process submitting jobs should be run
either as a batch job or on an interactive node - not on the biowulf
login node.










