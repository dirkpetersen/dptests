

document.querySelector('title').textContent = 'MitoSuite on Biowulf';
MitoSuite on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
 |



MitoSuite is a graphical tool for human mitochondrial genome profiling in massively parallel sequencing.



### References:


* Ishiya K, Ueda S. (2017) MitoSuite: a graphical tool for human mitochondrial genome profiling in massive parallel sequencing. PeerJ 5:e3406.


Documentation
* [MitoSuite Guide](https://mitosuite.com/#startguide)


Important Notes
* Module Name: mitosuite (see [the modules page](/apps/modules.html) for more information)
* Example files in /fdb/mitosuite



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.

Because MitoSuite is a graphical application, make sure to use X11-forwarding or NoMachine when [connecting to the login node](/docs/connect.html).

Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load mitosuite**
[user@cn3144 ~]$ **mitosuite**
(mitosuite window opens -- follow example in the guide at mitosuite.com)
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```







