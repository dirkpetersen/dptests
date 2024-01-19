

document.querySelector('title').textContent = 'Comsol on Biowulf';
Comsol on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
 |



Comsol Multiphysics (previously named Femlab) is a modeling package for the simulation of any
physical process you can describe with partial differential equations (PDEs). It features
state-of-the-art solvers that address complex problems quickly and accurately, while its 
intuitive structure is designed to provide ease of use and flexibility.



Documentation
* [Comsol: Support](http://www.comsol.com/support)


Important Notes
* Module Name: Comsol (see [the modules page](/apps/modules.html) for more information)
* Interactive Only



An [X-Windows](/docs/connect.html) connection is required.


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=14 --mem=20g --gres=lscratch:100**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3920 are ready for job

[user@cn3920 ~]$ **module load Comsol**
[user@cn3920 ~]$ **comsol -tmpdir /lscratch/$SLURM\_JOB\_ID**

[user@cn3920 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

You should see the Comsol Desktop image as below:


![Comsol Desktop](Comsol_53a.png)


After the first time running, the default location of Comsol output is generated in your home directory:



```
[user@biowulf]$ **ls ~/.comsol**
v61
```

Because your /home directory is limited in space, it is a good idea to move and symlink the ~/.comsol directory to your /data directory, where there is a much larger
amount of disk space available. This can be done like this:



```
[user@biowulf]$ **mv ~/.comsol /data/$USER/comsol ; ln -s /data/$USER/comsol ~/.comsol**
```





