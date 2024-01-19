

document.querySelector('title').textContent = 'mageck-vispr on Biowulf';
mageck-vispr on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



The mageck-vispr toolkit is used to QC, analyse, and visualize results from CRISPR/Cas9
screens. It includes the
* mageck algorithm for identifying essential genes from CRISPER/Cas9 screens
* vispr visualization platform for mageck results





The mageck workflow is implemented as a snakemake pipeline and runs automatically.
Vispr on the other hand is essentially a web application that will run a temporary server on
a compute node and the user will connect to it using a browser on his/her own computer
through an ssh tunnel.


### References:


* W. Li, H. Xu, T. Xiao, L. Cong, M. Love, F. Zhang, R. Irizarry, J. Liu, 
 M. Brown and S. Liu. *MAGeCK enables robust identification of essential
 genes from genome-scale CRISPR/Cas9 knockout screens*. Genome
 Biology 2014, 15:554.
 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/25476604) | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4290824/) | 
 [Journal](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0554-4)


Documentation
* [visper repo](https://bitbucket.org/liulab/vispr)
* [mageck repo](https://bitbucket.org/liulab/mageck-vispr/overview)


Important Notes
* Module Name: mageck-vispr (see [the modules page](/apps/modules.html) for more information)
* Example files in `$MAGECK_VISPR_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an interactive session with [sinteractive](/docs/userguide.html#int)
and use as shown below. See our [tunneling page](https://hpc.nih.gov/docs/tunneling/)
for more details on setting up ssh tunnels.



```

biowulf$ **sinteractive --cpus-per-task=8 --mem=16g --tunnel**
salloc.exe: Pending job allocation 31864544
salloc.exe: job 31864544 queued and waiting for resources
salloc.exe: job 31864544 has been allocated resources
salloc.exe: Granted job allocation 31864544
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn2692 are ready for job
srun: error: x11: no local DISPLAY defined, skipping

Created 1 generic SSH tunnel(s) from this compute node to
biowulf for your use at port numbers defined
in the $PORTn ($PORT1, ...) environment variables.


Please create a SSH tunnel from your workstation to these ports on biowulf.
On Linux/MacOS, open a terminal and run:

    ssh  -L 46763:localhost:46763 biowulf.nih.gov

For Windows instructions, see https://hpc.nih.gov/docs/tunneling

[user@cn3144]$ **module load mageck-vispr**

```

Copy test data and run mageck pipeline



```

[user@cn3144]$ **cp -r ${MAGECK\_VISPR\_TEST\_DATA}/esc-testdata .**
[user@cn3144]$ **cd esc-testdata**
[user@cn3144]$ **mageck-vispr init ./test\_workflow --reads \
 reads/ERR376998.subsample.fastq \
 reads/ERR376999.subsample.fastq \
 reads/ERR377000.subsample.fastq**
[user@cn3144]$ **tree test\_workflow**
  test_workflow/
  |-- [user   2.9K]  config.yaml
  |-- [user   7.3K]  README.txt
  `-- [user   5.8K]  Snakefile
  
  0 directories, 3 files

[user@cn3144]$ **cd test\_workflow**

```

Before running the workflow it is necessary to edit the automatically
generated config file. The generated file contains many comments. Here is
the edited file with the comments stripped for readability:



```

library: ../yusa_library.csv
species: mus_musculus
assembly: mm10

targets:
    genes: true
sgrnas:
    update-efficiency: false
    trim-5: AUTO
    len: AUTO

samples:
    esc1:
        - ../reads/ERR376999.subsample.fastq
    esc2:
        - ../reads/ERR377000.subsample.fastq
    plasmid:
        - ../reads/ERR376998.subsample.fastq

experiments:
    "ESC-MLE":
        designmatrix: ../designmatrix.txt

# the following setting is required for 0.5.6
correct_cnv: false

```

Once the config file has been modified to reflect the experimental
design, run the pipeline. Note that snakemake is used to run this locally, not
by submitting tasks as cluster jobs. Note that the snakemake installation
used for mageck-vispr has been renamed to `mageck-vispr-snakemake`
to avoid interfering with the general use snakemake.



```

[user@cn3144]$ **mageck-vispr-snakemake --cores=$SLURM\_CPUS\_PER\_TASK**

```

Next, start the vispr server for visualization. Note that sinteractive sets the
`$PORT1` variable to the port selected for tunneling.



```

[user@cn3144]$ **cd test\_workflow**
[user@cn3144]$ **vispr server --port $PORT1 --host=localhost results/ESC-MLE.vispr.yaml**
Loading data.
Starting server.

Open:  go to localhost:39335 in your browser.
Note: Safari and Internet Explorer are currently unsupported.
Close: hit Ctrl-C in this terminal.


```

On your local workstation, create an ssh tunnel to biowulf as describe in our
[tunneling documentation](https://hpc.nih.gov/docs/tunneling/).




![vispr web app](/images/mageck_vispr_screenshot.png)


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. mageck-vispr.sh), which uses the input file 'mageck-vispr.in'. For example:


Create a batch script for an existing config file similar to the following example:



```

#! /bin/bash
# this file is mageck.batch

module load mageck-vispr/0.5.4 || exit 1
cd /path/to/workdir

mageck-vispr-snakemake --cores=$SLURM_CPUS_PER_TASK

```

Submit to the queue with [sbatch](/docs/userguide.html):



```

biowulf$ **sbatch --cpus-per-task=8 --mem=16g mageck.batch**

```







