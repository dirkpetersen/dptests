

document.querySelector('title').textContent = 'RStudio Server on Biowulf';
RStudio Server on Biowulf
[![RStudio Logo](/images/RStudio-Ball.png)](http://www.rstudio.com)
[RStudio](https://www.rstudio.com/) is an integrated development
environment (IDE)for R. It includes a console, syntax-highlighting editor that
supports direct code execution, as well as tools for plotting, history,
debugging and workspace management.


This is RStudio Server, a browser-based
interface very similar to the standard RStudio desktop environment. RStudio Server
can be much more responsive and a generally-better experience when used remotely,
especially over a VPN.
 RStudio Server is
used interactively and requires an SSH tunnel connection
to biowulf. Please use an [sinteractive tunnel](/docs/tunneling) 
to create the tunnel necessary to access RStudio Server.



 This module is currently experimental and may change or disappear without warning. Please contact [staff@hpc.nih.gov](mailto:staff@hpc.nih.gov) if you run into issues while following the below steps. 



Common pitfalls
[top](#top)


**working directory is /home**
By default, the working directory is /home no matter where the interactive session is started. This could be changed by setwd(), which will set the workding directory to desired location. 

setwd("/data/somewhere")

  

Another way (more recommended) is using [project](https://support.posit.co/hc/en-us/articles/200526207-Using-RStudio-Projects), which will set it's own working directory.
 File -> New Project -> New Directory -> Create New Project (assign directory name, and choose the directory location such as "/data/XXX/")
  



**My R sessions are/aren't being restored when I reconnect**
R sessions inside RStudio server are persistent within a single sinteractive session.
 Different sinteractive jobs cannot share R sessions and job sessions are not persistent
 between jobs. Please save your work to disk or explicitly save your session to an RData
 file by manually saving the workspace with the Session -> Save Workspace As... menu.
 


**I am getting a permission denied error**
R sessions inside RStudio server are persistent within a single sinteractive session.
 Every time you run rstudio-server, you will receive a new randomized password. Please make
 sure to use the new link or password each time you restart the server.
 



Running RStudio Server interactively
Allocate an [sinteractive](https://hpc.nih.gov/docs/userguide.html#int)
session. It is required to allocate at least a small amount of [lscratch](/docs/userguide.html#local) for temporary storage for R.


It is also required to set up an [sinteractive tunnel](/docs/tunneling).




```

[user@biowulf ~]$ **sinteractive --mem=10g --gres=lscratch:5 --tunnel**
salloc.exe: Pending job allocation 15323416
salloc.exe: job 15323416 queued and waiting for resources
salloc.exe: job 15323416 has been allocated resources
salloc.exe: Granted job allocation 15323416salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn1640 are ready for job

Created 1 generic SSH tunnel(s) from this compute node to
biowulf for your use at port numbers defined
in the $PORTn ($PORT1, ...) environment variables.


Please create a SSH tunnel from your workstation to these ports on biowulf.
On Linux/MacOS, open a terminal and run:

    ssh  -L 39689:localhost:39689 user@biowulf.nih.gov

For Windows instructions, see https://hpc.nih.gov/docs/tunneling
[user@cn1640 ~]$ **module load rstudio-server**
[+] Loading gcc  9.2.0  ...
[-] Unloading gcc  9.2.0  ...
[+] Loading gcc  9.2.0  ...
[+] Loading openmpi 4.0.5  for GCC 9.2.0
[+] Loading ImageMagick  7.0.8  on cn4280
[+] Loading HDF5  1.10.4
[-] Unloading gcc  9.2.0  ...
[+] Loading gcc  9.2.0  ...
[+] Loading NetCDF 4.7.4_gcc9.2.0
[+] Loading pandoc  2.17.1.1  on cn4280
[+] Loading pcre2 10.21  ...
[+] Loading R 4.2.2
[+] Loading rstudio-server  2023.03.0-386
[user@cn1640 dir]$ **rstudio-server**

Please ensure you have set up the SSH port forwarding as described in the sinteractive instructions.

Please connect to http://localhost:39689/auth-sign-in?user=test2&password=nRmzfPWh_X8Z-03hbDjPz3bm
Use your username 'user' and the pasword 'nRmzfPWh_X8Z-03hbDjPz3bm' to login


```

If you have set up the SSH tunnel with the command or instructions provided (please make sure to
use the specific port/command given by *your* interactive session) you should be able to
load the link and have RStudio Server load in your browser.


If you are presented with a login screen, please use your username and the password provided


Documentation
* [RStudio Web Site](http://www.rstudio.com)* [RStudio Resources](https://www.rstudio.com/online-learning/#R)* [RStudio Server Web Site](https://posit.co/products/open-source/rstudio-server/)






