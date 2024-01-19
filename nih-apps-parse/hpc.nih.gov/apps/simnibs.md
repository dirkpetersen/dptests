

document.querySelector('title').textContent = 'SimNIBS on Biowulf';
SimNIBS on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Quick start](#quick)
 |



SimNIBS is used to simulate non-invasive brain stimulation. Calculations of the 
electric field induced by transcranial magnetic stimulation (TMS) and transcranial
direct current stimulation (tDCS) are supported.



### References:


* A. Thielscher, A. Antunes, and G. B. Saturnino. *Field modeling for
 transcranial magnetic stimulation: a useful tool to understand the
 physiological effects of TMS?*. Conf Proc IEEE Eng Med Biol Soc, 2015.
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/26736240) | 
 [Journal](http://ieeexplore.ieee.org/xpls/icp.jsp?arnumber=7318340)


Documentation
* [Home page](https://simnibs.github.io/simnibs/build/html/index.html)


Important Notes
* Module Name: simnibs (see [the modules page](/apps/modules.html) for more information)
* Example files in `$SIMNIBS_TEST_DATA`



Quickstart
Copy the test data and set up your environment



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load simnibs**
[+] Loading simnibs 3.0 on cn3144 

[user@cn3144]$ **cd /data/$USER**
[user@cn3144]$ **cp $SIMNIBS\_TEST\_DATA/simnibs-3.0-examples.zip .**
[user@cn3144]$ **unzip simnibs-3.0-examples.zip**

```

The first step would usually be the construction of a 3D mesh from MRI
images. Since this step is time consuming, it should be run as a batch job. You have two options: 
**mri2mesh** and **headreco**. **mri2mesh** uses FSL for skull segmentation and FreeSurfer for 
brain segmentation. You would create a batch script similar to the following



```

#!/bin/bash
#SBATCH --job-name=mri2mesh-job
#SBATCH --cpus-per-task=8
#SBATCH --mem=10g
#SBATCH --time=24:00:00

module load simnibs fsl freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh

cd /data/$USER/simnibs_examples
mkdir new_mesh
cd new_mesh
mri2mesh --all ernie ../ernie/org/ernie_T1.nii.gz ../ernie/org/ernie_T2.nii.gz

```

and submit it as a batch job with



```

[user@biowulf]$ **sbatch make\_mesh\_mri2mesh.sh**

```

If you wish to use **headreco** (which uses matlab, SPM, and CAT12), you can do so with a 
batch script as follows:



```

#!/bin/bash
#SBATCH --job-name=headreco-job
#SBATCH --cpus-per-task=32
#SBATCH --mem=8g
#SBATCH --time=7:00:00

export JAVA_OPTS="-Xmx100g -XX:-UseGCOverheadLimit"
export _JAVA_OPTIONS="-Xmx100g -XX:-UseGCOverheadLimit"
export JAVA_TOOL_OPTIONS="-Xmx100g -XX:-UseGCOverheadLimit"

cd /data/$USER/simnibs_examples/ernie

module load simnibs/3.0 matlab/2018b

headreco all ernie org/ernie_T1.nii.gz org/ernie_T2.nii.gz

```

and submit it as a batch job with



```

[user@biowulf]$ **sbatch make\_mesh\_headreco.sh**

```

Since the example directory already contains a mesh reconstruction of
the example data, we can proceed with the tutorial without waiting
for the batch job to finish. This next part uses the SimNIBS GUI and
we will run it on an interactive node. Since this is a GUI program,
using [Nomachine NX](https://hpc.nih.gov/docs/nx.html) to connect 
to biowulf is likely to give better resonsiveness.



```

[user@biowulf]$ **sinteractive --mem=10g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load simnibs**
[+] Loading simnibs  3.0  on cn3144 

[user@biowulf]$ **cd /data/$USER/simnibs\_test**
[user@biowulf]$ **simnibs\_gui &**

```

This will start the following GUI:


![GUI 1](/images/simnibs_gui.png)
 You will find a step-by-step guide on how to set up tDCS and TMS simulations (and more) at the [Simnibs GUI Tutorial](https://simnibs.github.io/simnibs/build/html/tutorial/gui.html).






