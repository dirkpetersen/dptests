

document.querySelector('title').textContent = 'Abyss on HPC';
Abyss on HPC


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
  |

  ABySS is a de novo, parallel, paired-end sequence assembler that is designed 
 for short reads. The parallel version is implemented using MPI and is capable 
 of assembling larger genomes. 
 



### References:

 * <https://genome.cshlp.org/content/27/5/768>


Documentation * <https://github.com/bcgsc/abyss>



Important Notes * Module Name: abyss (see [the modules 
 page](/apps/modules.html) for more information)
* MPI app; when running abyss on multiple nodes, please follow the syntax in the example below instead of using traditional MPI syntax.





Batch job (in MPI mold)
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. batch.sh). For example:



```

#!/bin/bash
#SBATCH --job-name="abyss"
#SBATCH --mail-type=BEGIN,END
### sbatch --partition=multinode --ntasks=16 --nodes=2 --time=24:00:00 --mem=60g batch.sh

cd /data/$USER/abyss
module load abyss
`which abyss-pe` np=${SLURM_NTASKS} j=8 k=25 n=10 in='/data/$USER/File_1.fq /data/$USER/File_2.fq' name=OutputPrefix

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --partition=multinode --ntasks=16 --nodes=2 --time=24:00:00 --mem=60g batch.sh
```

 The job runs 16 tasks ( np=${SLURM\_NTASKS} ) on 2 nodes with 8 (j=8) cpus each node.

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1;abyss-pe np=${SLURM_NTASKS} j=8 k=25 n=10 in='1.fq 2.fq' name=out
cd dir2;abyss-pe np=${SLURM_NTASKS} j=8 k=25 n=10 in='1.fq 2.fq' name=out
cd dir3;abyss-pe np=${SLURM_NTASKS} j=8 k=25 n=10 in='1.fq 2.fq' name=out

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm --module abyss --sbatch "--partition=multinode --ntasks=16 --nodes=2 --time=24:00:00 --mem=60g"
```

 where 
 

|  |  |
| --- | --- |
| --sbatch | use this flag to pass sbatch flags to swarm |
| --module  | Loads the module for each subjob in the swarm  |







