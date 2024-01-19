

document.querySelector('title').textContent = 'MIPAV on Biowulf';
MIPAV on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
 |



[![](http://mipav.cit.nih.gov/about.asp_files/splash.gif)](http://mipav.cit.nih.gov)



The MIPAV (Medical Image Processing, Analysis, and Visualization) application enables quantitative analysis and visualization of medical images of 
numerous modalities such as PET, MRI, CT, or microscopy. Using MIPAV's standard user-interface and analysis tools, researchers at remote sites 
(via the internet) can easily share research data and analyses, thereby enhancing their ability to research, diagnose, monitor, and treat medical disorders. 




Documentation
* [MIPAV website](http://mipav.cit.nih.gov)


Important Notes
* Module Name: MIPAV (see [the modules page](/apps/modules.html) for more information)
* Don't run MIPAV on the login node. Use an interactive session instead as in the example below.
* Singlethreaded/multithreaded.
* Sample data in /usr/local/apps/MIPAV/sample\_data/* Environment variables set 
	+ MIPAV\_HOME



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

[user@cn3144 ~]$ **module load MIPAV**

[user@cn3144 ~]$  **mipav &**
![](/images/mipav1.jpg)

![](/images/mipav2.jpg)

![](/images/mipav3.png)

![](/images/mipav4.png)
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$









```










