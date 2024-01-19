

document.querySelector('title').textContent = 'DYNAMO on Biowulf';
DYNAMO on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
 |



Dynamo is a software environment for subtomogram averaging of cryo-EM data. 



### References:


* [Dynamo: A flexible, user-friendly development tool for subtomogram averaging of cryo-EM data in high-performance computing environments](https://www.sciencedirect.com/science/article/pii/S1047847711003650). Castaño-Díez D, Kudryashev M, Arheit M, Stahlberg H., J Struct Biol. 2012.


Documentation
* [Dynamo Main Site](https://wiki.dynamo.biozentrum.unibas.ch/w/index.php/Main_Page)


Important Notes
* Module Name: dynamo (see [the modules page](/apps/modules.html) for more information)



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --gres=gpu:k80:2,lscratch:10 --mem=20g -c14**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load dynamo**
[+] Loading dynamo  1.4.01  on cn3144 

[user@cn3144 ~]$ **source dynamo\_activate\_linux\_shipped\_MCR.sh /usr/local/apps/dynamo/1.4.01/**
MCR for linux has been  found in location: /dynamo/1.4.01/MCRLinux/v94
Activating Dynamo as standalone
-------------------------------------------------------------
 
              Dynamo Setup for standalone Linux 
 
-------------------------------------------------------------
 
       
1- Adding execution paths
   ----------------------
       
    Identifying directories of Dynamo distribution hanging from the specified Dynamo root directory: 
     /dynamo/1.4.01
 
    [OK] Done updating Linux paths.
       
2- Testing OpenMP
   --------------
       
    checking functionality of OpenMP wrappers....

    [OK] OpenMP tools seem functional
       
3- Testing MCR libraries
   ---------------------
       
    [OK] MCR library seems active
  
  
------------------------------------
Dynamo is ready for standalone use in your Linux system
 
[user@cn3144 ~]$ **export MCR\_CACHE\_ROOT=/lscratch/$SLURM\_JOBID**

[user@cn3144 ~]$ **dynamo**
Initializing MATLAB Runtime version 9.4
"MCR libraries starting dynamo. Starting the Dynamo console can take some time. A good choice of MCR_CACHE_ROOT might speed this up."
Launching the console....
 (use "dynamo x" for a Dynamo graphics terminal  on a independent window)
 
     Dynamo console for standalone modus
     -----------------------------------
       
       
     * Commands starting with the symbol "!" will be passed to the system
     * Commands starting with the symbol "\" will be executed as Matlab expressions
     * Strings including the symbol ":" are treated as database queries
     * Graphic output: enabled
     * Type "help" for a list of console control commands
       
       
Dynamo > **exit**
ok, bye!
 

Exiting Dynamo console.
Returning control to system.
 

 
 Note after exiting the Dynamo console:
   
       if you experience difficulties with the text display you might need to type: 
       stty echo 
       in your system shell to fix them.
 
[delete:pmodelpool.Manager] Deletion of model pool manager from memory.

[user@cn3144 ~]**dynamo x &**
[1] 5233
[user@cn3144 ~]$ Initializing MATLAB Runtime version 9.4
"MCR libraries starting dynamo. Starting the Dynamo console can take some time. A good choice of MCR_CACHE_ROOT might speed this up."
Launching a Dynamo graphics terminal in a separate window....


```

The above should display dynamo's graphical interface.







