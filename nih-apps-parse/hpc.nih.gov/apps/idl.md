

document.querySelector('title').textContent = 'IDL/ENVI on Biowulf';
IDL/ENVI on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
 |



IDL provides powerful, core visualization and analysis functionality, and capabilities that allow data analysts and developers to leverage IDL's power in multiple software environments.


ENVI image analysis software is used to extract meaningful information from imagery to make better decisions. 



### References:


* Paper


Documentation
* [IDL and ENVI documentation](https://www.harrisgeospatial.com/docs/)


Important Notes
* Module Name: IDL or ENVI(see [the modules page](/apps/modules.html) for more information)
* IDL and ENVI cannot be run on the Biowulf login node. Instead, allocate an interactive session as described, and run IDL/ENVI there.
* Licenses: There are 20 IDL licenses on Biowulf. Any single job requires 6 licenses, so at most 3 sessions can be run simultaneously with the Biowulf IDL licenses. 
* Example files available in /usr/local/apps/Exelis/8.5\_5.3/idl/examples



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):


To run the IDL GUI, you should have connected to Biowulf with an Xwindows session. (see the [Connect to Biowulf](/docs/connect.html) page for more info)

```

[user@biowulf]$ **sinteractive --license=idl:6**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load IDL**

[user@cn3144 ~]$  **idl**    #command-line session
IDL> x = FINDGEN(41)/10 - 2
IDL> gauss = EXP(-x^2)
IDL> myPlot = BARPLOT(x, gauss, TITLE='Gaussian Distribution', XTITLE='$\it x$', YTITLE='$\it f(x)$', YRANGE=[0,1.1])
IDL> myText = TEXT(0.75,0.85,'$\it f(x)=exp(-x^2)$', /DATA, FONT_SIZE=24)
IDL>myPlot.Save, "gaussian.png", BORDER=10, RESOLUTION=300, /TRANSPARENT
IDL> exit

![](/images/idl_gaussian.png)

[user@cn3144 ~]$ **idlde** # GUI interface

![](/images/idl1.jpg)  

![](/images/idl2.jpg)


[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```










