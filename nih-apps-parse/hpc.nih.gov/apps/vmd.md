

document.querySelector('title').textContent = 'VMD';
VMD


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive VMD on Biowulf](#int)
[Plugins](#plugins)
[On your desktop using hpcdrive](#hpcdrive) 
 |


![VMD logo](/images/VMD_logo_left.gif)
VMD is a molecular visualization program for displaying, 
 animating, and analyzing large biomolecular systems using 3-D 
 graphics and built-in scripting. It has powerful and 
 comprehensive filtering and configuration capabilities. It is
 especially well-suited for analyzing [NAMD](/apps/NAMD.html) results.




Documentation
* [VMD Documentation](http://www.ks.uiuc.edu/Research/vmd/current/docs.html)




Important Notes
* Module Name: vmd (see [the modules page](/apps/modules.html) for more information)
* Number of CPUs: By default, VMD will try to use all CPUs on the system. The VMD module on Biowulf has been 
written to set the number of CPUs to the number of allocated CPUs using the environment variable VMDFORCECPUCOUNT.




Running VMD in an interactive session on Biowulf

 VMD is a graphics program, and therefore you need a graphics connection to Biowulf. We recommend
 installing NX ([Installation Instructions](https://hpc.nih.gov/docs/nx.html)). 
 
 Once you have NX installed, open a session to Biowulf. Confirm that the graphics connection is working by typing 'xclock' -- you should
 see a clock appear. 
 
**Please do not run VMD on the login node,** which is shared by many users.
 Instead, start an [interactive session](/docs/userguide.html#int) and run the program. 
 If you expect to be running calculations with VMD, you may want to request more than the default 2 CPUs. In that case,
 use the flag --cpus-per-task=# in your sinteractive command.
 
 Sample VMD session: (user input in bold)
 
```

biowulf% **sinteractive --cpus-per-task=4** 
	
salloc.exe: Pending job allocation 13212187
salloc.exe: job 13212187 queued and waiting for resources
salloc.exe: job 13212187 has been allocated resources
salloc.exe: Granted job allocation 13212187
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0851 are ready for job

[user@cn0851 ~]$ **module load vmd**
[+] Loading netpbm 10.86 on cn0851 
[+] Loading vmd v 1.9.3  ... 
       Using 4 CPUs

[user@cn0851 ~]$ **vmd**
/usr/local/apps/vmd/1.9.3/lib/vmd/vmd_LINUXAMD64: /lib64/libGL.so.1: no version information available (required by /usr/local/apps/vmd/1.9.3/lib/vmd/vmd_LINUXAMD64)
Info) VMD for LINUXAMD64, version 1.9.3 (November 30, 2016)
Info) http://www.ks.uiuc.edu/Research/vmd/                         
Info) Email questions and bug reports to vmd@ks.uiuc.edu           
Info) Please include this reference in published work using VMD:   
Info)    Humphrey, W., Dalke, A. and Schulten, K., `VMD - Visual   
Info)    Molecular Dynamics', J. Molec. Graphics 1996, 14.1, 33-38.
Info) -------------------------------------------------------------
Info) Multithreading available, 4 CPUs detected.
Info)   CPU features: SSE2 AVX AVX2 FMA KNL:AVX-512F+CD+ER+PF 
Info) Free system memory: 271GB (71%)
Info) No CUDA accelerator devices available.
Warning) Detected X11 'Composite' extension: if incorrect display occurs
Warning) try disabling this X server option.  Most OpenGL drivers
Warning) disable stereoscopic display when 'Composite' is enabled.
Info) OpenGL renderer: llvmpipe (LLVM 5.0, 256 bits)
Info)   Features: STENCIL MDE CVA MTX NPOT PP PS GLSL(OVF) 
Info)   Full GLSL rendering mode is available.
Info)   Textures: 2-D (8192x8192), 3-D (512x512x512), Multitexture (8)
Info) Dynamically loaded 2 plugins in directory:
Info) /usr/local/apps/vmd/1.9.3/lib/vmd/plugins/LINUXAMD64/molfile

```


At this point you should see two VMD windows appear. VMD is menu-driven, so now you can use the VMD File menu to load a PDB or traajectory file.

![](/images/vmd1.png)

![](/images/vmd2.png)

![](/images/vmd3.png)







Plugins
There are two plugins available for VMD, [Molcontroller](https://github.com/ChenchenWu-hub/Molcontroller) and [RIP-MD](https://github.com/DLab/RIP-MD).


To use these plugins, you will need to add the following lines to your ~/.vmdrc file:



```
lappend auto_path /usr/local/apps/vmd/Molcontroller/molcontrol1.0/
vmd_install_extension molcontrol molcontroller_tk "Modeling/Molcontrol"

lappend auto_path /usr/local/apps/vmd/RIP-MD/
vmd_install_extension ripmd ripmd_tk "Analysis/RIP-MD"
```

The Molcontroller plugin appears in the VMD Main panel under *Extensions --> Modeling --> Molcontrol*.



RIP-MD appears in the VMD Main panel under *Extensions --> Analysis --> RIP-MD*.





For additional plugins, please contact *staff@hpc.nih.gov*.



Running VMD on your desktop using files on Biowulf

1. Map your Biowulf /home or /data area onto your desktop computer using [hpcdrive](/docs/hpcdrive.html). 
- Install VMD on your desktop computer. [[VMD Downloads](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD)]
 [[VMD Installation Guide](http://www.ks.uiuc.edu/Research/vmd/current/ig/ig.html)]
- Now you can start up VMD on your desktop, then load a file from your Biowulf /data or /home area.


















