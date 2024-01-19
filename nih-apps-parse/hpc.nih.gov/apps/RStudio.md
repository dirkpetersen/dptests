

document.querySelector('title').textContent = 'RStudio on Biowulf';
RStudio on Biowulf
[![RStudio Logo](/images/RStudio-Ball.png)](http://www.rstudio.com)
[RStudio](https://www.rstudio.com/) is an integrated development
environment (IDE)for R. It includes a console, syntax-highlighting editor that
supports direct code execution, as well as tools for plotting, history,
debugging and workspace management. RStudio is 
used interactively and requires a graphical connection
to biowulf. Please use the [NoMachine NX client](/docs/nx.html) 
to create the graphical connection as it is the most reliable for RStudio.



 There is RStudio Server (a browser-based interface) available with better performance:
 [rstudio-server](https://hpc.nih.gov/apps/rstudio-server.html)



Common pitfalls
[top](#top)


**working directory is /home**
By default, the working directory is /home no matter where the interactive session is started. This could be changed by setwd(), which will set the workding directory to desired location. 

setwd("/data/somewhere")

  

Another way (more recommended) is using [project](https://support.posit.co/hc/en-us/articles/200526207-Using-RStudio-Projects), which will set it's own working directory.
 File -> New Project -> New Directory -> Create New Project (assign diretory name, and choose the directory location such as "/data/XXX/")
  


  

**XDG\_RUNTIME\_DIR not set**
This is a harmless warning:  
 QStandardPaths: XDG\_RUNTIME\_DIR not set, defaulting to '/XXXX/'

  

**Can't open Rstudio** 
Rstudio requires graphical connection, please use [NX](/docs/nx.html) to avoid these errors:
 
 libGL error: No matching fbConfigs or visuals found
   
libGL error: failed to load driver: swrast
   
Unrecognized OpenGL version
 or 
 
 qt.qpa.xcb: could not connect to display 
   
qt.qpa.plugin: Could not load the Qt platform plugin "xcb" in "" even though it was found.
 



Running RStudio interactively
Connect to Biowulf with NoMachine ([NX](/docs/nx.html)), then open a terminal inside of NX.


![rstudio window](/images/rstudio_nx1.png)


![rstudio window](/images/rstudio_nx2.png)

Allocate an [sinteractive](https://hpc.nih.gov/docs/userguide.html#int)
session. It is generally recommended to allocate at least a small amount of [lscratch](/docs/userguide.html#local) for temporary storage for R.




```

[user@biowulf ~]$ **sinteractive --mem=10g --gres=lscratch:5**
salloc.exe: Pending job allocation 15323416
salloc.exe: job 15323416 queued and waiting for resources
salloc.exe: job 15323416 has been allocated resources
salloc.exe: Granted job allocation 15323416salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn1640 are ready for job
[user@cn1640 ~]$ **module load rstudio R**
[+] Loading Rstudio 1.1.447 
Remember to load an R module before starting Rstudio
...
[user@cn1640 dir]$ **rstudio &**

```


![rstudio window](/images/rstudio_nx3.png)

Documentation
* [RStudio Web Site](http://www.rstudio.com)* [RStudio Resources](https://www.rstudio.com/online-learning/#R)








