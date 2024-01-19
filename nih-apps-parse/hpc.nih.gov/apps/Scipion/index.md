

document.querySelector('title').textContent = 'Scipion';
Scipion


|  |
| --- |
| 
Quick Links
[Interactive job on Biowulf](#int)
[Batch job on Biowulf](#sbatch)
[Documentation](#doc)
 |



[Scipion](https://github.com/I2PC/scipion) is an image processing framework to obtain 3D models of macromolecular complexes using Electron Microscopy (3DEM). It integrates several software packages and presents an unified interface for both biologists and developers. Scipion allows to execute workflows combining different software tools, while taking care of formats and conversions. Additionally, all steps are tracked and can be reproduced later on.



There may be multiple versions of Scipion available. An easy way of selecting the version is to use [modules](/apps/modules.html). To see the modules available, type



```
module avail scipion
```

To select a module, type



```
module load scipion/[ver]
```

where [ver] is the version of choice.


### Environment variables set:


* $PATH


Configuration
Scipion is a package management framework. Much of what Scipion does is controlled via the local scipion.conf file. This is located in ~/.config/scipion:



```
module load scipion
scipion config
cat ~/.config/scipion/scipion.conf
```

This file is created when Scipion is first used. The **SCIPION\_USER\_DATA** parameter in the **DIRS\_LOCAL** section of this file should be edited to direct where generated files are written, as by default it will write to your /home directory:



```
[DIRS_LOCAL]
SCIPION_USER_DATA = /data/**$USER**/apps/scipion/ScipionUserData
SCIPION_LOGS = %(SCIPION_USER_DATA)s/logs
SCIPION_TMP = %(SCIPION_USER_DATA)s/tmp
```

**$USER** should be changed to your username.


Interactive job on Biowulf
Once an interactive session is allocated, e.g.



```
sinteractive --cpus-per-task=16 --mem=24g --gres=lscratch:50 --time=8:00:00
```

Scipion can be run like so:



```

module load scipion
scipion tests em.workflows.test_workflow_spiderMDA

```

Scipion can be run as a GUI. This requires an [X11 session](https://hpc.nih.gov/docs/connect.html).



```
scipion
```

![GUI](scipion_1.png)
A batch job can be launched from the GUI like this:


Double click on the project step:


![GUI Batch 1](scipion_gui_batch_1.png)
Making sure all the input parameters are correct, change 'Use queue?' from 'No' to 'Yes', then click 'Execute':


![GUI Batch 2](scipion_gui_batch_2.png)
The Queue parameters window will appear. Choose the appropriate partition (queue):


![GUI Batch 3](scipion_gui_batch_3.png)
Make sure the memory, time and local scratch space are adequate, then click 'OK':


![GUI Batch 4](scipion_gui_batch_4.png)
Batch job on Biowulf
Create a batch input file (e.g. scipion.sh):



```

#!/bin/bash
module load scipion
scipion tests em.workflows.test_workflow_spiderMDA

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=16 --mem=24g scipion.sh
```

Documentation
* Scipion Main Site: <https://github.com/I2PC/scipion>




