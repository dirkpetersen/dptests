

document.querySelector('title').textContent = 'MaxQuant on Biowulf';
MaxQuant on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
 |



MaxQuant is a quantitative proteomics software package designed for analyzing large mass-spectrometric data sets. It is specifically aimed at high-resolution MS data. Several labeling techniques as well as label-free quantification are supported.



### References:


* Cox J, Mann M.
 [MaxQuant enables high peptide identification rates, individualized p.p.b.-range mass accuracies and proteome-wide protein quantification.](https://www.ncbi.nlm.nih.gov/pubmed/19029910) 
*Nat Biotechnol. 2008 Dec;26(12):1367-72.** Tyanova S, Temu T, Cox J.
 [The MaxQuant computational platform for mass spectrometry-based shotgun proteomics.](https://www.ncbi.nlm.nih.gov/pubmed/27809316) 
*Nat Protoc. 2016 Dec;11(12):2301-2319.*


Documentation
* [MaxQuant Main Site](https://www.maxquant.org/)


Important Notes
* Module Name: maxquant (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded
 * This application produces HTML reports. You can use [hpcdrive](/docs/hpcdrive.html) to view these reports on your local workstation.


This application requires an [X-Windows connection](/docs/connect.html). Users are encouraged to use [NX](https://hpc.nih.gov/docs/nx.html) as their X11 servers.


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

[user@cn3144 ~]$ **module load maxquant**
[user@cn3144 ~]$ **MaxQuantGui.exe**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

An interactive window will pop up:


![MaxQuant window](maxquant.PNG)


