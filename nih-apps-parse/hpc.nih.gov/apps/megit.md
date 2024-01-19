

document.querySelector('title').textContent = "TEMPLATE";
MeGit on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
 |



Egit is a plug-in for the Eclipse IDE to manage git repositories. Minimal Egit (MeGit) is a standalone version of Egit. MeGit is a graphical front-end for handling git repositories, git branch operations, analyzing git history, etc.



Documentation
* [MeGit Main Site](https://github.com/eclipsesource/megit)
* [EGit User Guide](https://wiki.eclipse.org/EGit/User_Guide)


Important Notes
This application requires a [graphical connection using NX](/docs/connect.html#nx)


* Module Name: megit (see [the modules page](/apps/modules.html) for more information)
 * This module is available in interactive jobs and on Helix.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
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

[user@cn3144 ~]$ **module load megit**
[+] Loading java 17.0.3.1  ...
[+] Loading megit 0.4.0

[user@cn3144 ~]$ **megit**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226

[user@biowulf ~]$

```








