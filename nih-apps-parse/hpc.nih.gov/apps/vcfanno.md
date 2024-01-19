

document.querySelector('title').textContent = 'Vcfanno on HPC';
Vcfanno on HPC
vcfanno annotates a VCF with any number of *sorted* and tabixed 
 input BED, BAM, and VCF files in parallel. It does this by finding overlaps 
 as it streams over the data and applying user-defined operations on the 
 overlapping annotations. 


In order to parallelize, work is broken down as follows. A slice (array) 
 of query intervals is accumulated until a specified number is reached (usually 
 ~5K-25K) or a gap cutoff is exceeded; at that point, the bounds of the region 
 are used to perform a tabix (or any regional) query on the database files. 
 This is all done in [irelate](https://github.com/brentp/irelate). 
 vcfanno then iterates over the streams that result from the tabix queries 
 and finds intersections with the query stream. This is a parallel chrom-sweep. 
 This method avoids problems with chromosome order.


For VCF, values are pulled by name from the INFO field. For BED, values 
 are pulled from (1-based) column number. For BAM, depth (count), "mapq" 
 and "seq" are currently supported.


### 


Documentation
* <https://github.com/brentp/vcfanno>



Important Notes
* Module Name: vcfanno (see [the 
 modules page](/apps/modules.html) for more information)
* Example files: /usr/local/apps/vcfanno/version/example/





Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load vcfanno**
[user@cn3144 ~]$ **vcfanno -lua example/custom.lua example/conf.toml example/query.vcf.gz**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. batch.sh). For example:



```

#!/bin/bash
set -e
module load vcfanno
vcfanno -lua example/custom.lua example/conf.toml example/query.vcf.gz
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; vcfanno -lua example/custom.lua example/conf.toml example/query.vcf.gz
cd dir2; vcfanno -lua example/custom.lua example/conf.toml example/query.vcf.gz
cd dir3; vcfanno -lua example/custom.lua example/conf.toml example/query.vcf.gz

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm --module vcfanno
```



