

document.querySelector('title').textContent = 'pdf2svg on Biowulf';
pdf2svg on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



pdf2svg is a simple tool for converting PDF files to SVG format.


Documentation
* [pdf2svg Main Site](https://github.com/dawbarton/pdf2svg)


Important Notes
* Module Name: `pdf2svg` (see [the modules page](/apps/modules.html) for more information)



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session (user input in **bold**):



```

[user@biowulf ~]$ **sinteractive --gres=lscratch:10**
salloc.exe: Pending job allocation 11085118
salloc.exe: job 11085118 queued and waiting for resources
salloc.exe: job 11085118 has been allocated resources
salloc.exe: Granted job allocation 11085118
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0848 are ready for job
srun: error: x11: no local DISPLAY defined, skipping
error: unable to open file /tmp/slurm-spank-x11.11085118.0
slurmstepd: error: x11: unable to read DISPLAY value

[user@cn0848 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn0848 11085118]$ **wget https://upload.wikimedia.org/wikipedia/commons/d/d3/Test.pdf**

[user@cn0848 11085118]$ **module load pdf2svg**
[+] Loading pdf2svg  0.2.3  on cn0848

[user@cn0848 11085118]$ **pdf2svg Test.pdf test.svg**

[user@cn0848 11085118]$ **exit**
exit
salloc.exe: Relinquishing job allocation 11085118

[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. convert\_pdfs.sh). For example:



```

#!/bin/bash
set -e
module load pdf2svg
for pdf if ~/pdfs; do pdf2svg $pdf ~/svgs`basename -s .pdf $pdf`.svg; done

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch convert_pdfs.sh
```









