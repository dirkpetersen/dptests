

document.querySelector('title').textContent = 'cell2location on Biowulf';
cell2location on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
 |


 Cell2location: Comprehensive mapping of tissue cell architecture via integrated single cell and spatial transcriptomics (cell2location model)



### References:


* [Cell2location maps fine-grained cell types in spatial transcriptomics](https://pubmed.ncbi.nlm.nih.gov/35027729/). Nat Biotechnol. 2022 May;40(5):661-671. doi: 10.1038/s41587-021-01139-4. Epub 2022 Jan 13


Documentation
* <https://github.com/BayraktarLab/cell2location>
* <https://cell2location.readthedocs.io/en/latest/>


Important Notes
* Module Name: cell2location (see [the modules page](/apps/modules.html) for more information)
* */usr/local/apps/cell2location/test.ipynb* is available for demonstration which is from <https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html>.
 * This installation of cell2location is able to make use of gpu.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **module load tmux**
[user@biowulf]$ **tmux**
[user@biowulf]$ **sinteractive --time=8:00:00 --gres=lscratch:5,gpu:p100:1 --mem=5g --tunnel**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job
                                                                           
Created 1 generic SSH tunnel(s) from this compute node to                  
biowulf for your use at port numbers defined                               
in the $PORTn ($PORT1, ...) environment variables.                         
                                                                           
                                                                           
Please create a SSH tunnel from your workstation to these ports on biowulf.
On Linux/MacOS, open a terminal and run:                                   
                                                                           
    ssh  -L 33327:localhost:33327 biowulf.nih.gov                          
                                                                           
For Windows instructions, see https://hpc.nih.gov/docs/tunneling       


[user@cn4224 ~]$ **cd /lscratch/$SLURM\_JOBID**
[user@cn4224 46116226]$ **cp /usr/local/apps/cell2location/test.ipynb**

```


Now open a local terminal from your desktop (powershell for windows or terminal for Mac/Linux) and run the command below. Replace 33327 with the port number you see above.

```
ssh  -L 33327:localhost:33327 biowulf.nih.gov
```

Back to the interactive session, 


```

[user@cn4224 46116226]$ **module load cell2location**
[+] Loading cell2location 0.1.3 
[user@cn4224 46116226]$ **jupyter notebook --ip localhost --port $PORT1 --no-browser**

```


Copy and paste the URL from the output to your choice of internet browser and start jupyter notebook.  

Open the *test.ipynb* file and start the test.

```

[user@cn4224 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226

[user@biowulf ~]$

```







