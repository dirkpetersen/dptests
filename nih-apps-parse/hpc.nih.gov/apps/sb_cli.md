

document.querySelector('title').textContent = 'Seven Bridges CLI on Biowulf';
Seven Bridges CLI on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
 |



The Seven Bridges CLI, or SB CLI, lets users access a Seven Bridges platform using the API from the command line. The interface allows users to configure an API endpoint along with a corresponding authentication token. Use these credentials users can access a variety of information like their projects, tasks, and files, which can be downloaded or uploaded. 



Documentation
* [Seven Bridges CLI Command Reference](https://docs.sevenbridges.com/docs/general-commands)


Important Notes
* Module Name: sb\_cli (see [the modules page](/apps/modules.html) for more information)
* The module can be loaded on Helix or on compute nodes with, for example, an sinteractive session



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. The module is also available on Helix.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load sb\_cli**
[+] Loading sb_cli 0.18.4  ...

[user@cn3144 ~]$ **sb configure**
Seven Bridges API endpoint: https://cgc-api.sbgenomics.com
Seven Bridges authorization token: ****************
Seven Bridges advance access (true/false): false

[user@cn3144 ~]$ **sb --output table files list \
 --project username/test-project \
 --origin 1c46ecb3-772d-46b4-88db-6b423e75e811**
+--------------------------+-------------------------------+
|            ID            |             Name              |
+--------------------------+-------------------------------+
| 61d5f311e722267a17c521f6 | G20502.22Rv1.2_fastqc.b64html |
+--------------------------+-------------------------------+
| 61d5f2f9e722267a17c521b3 | G20502.22Rv1.2_fastqc.zip     |
+--------------------------+-------------------------------+

[user@cn3144 ~]$ **sb download --file 61d5f311e722267a17c521f6**
Downloading file `G20502.22Rv1.2_fastqc.b64html` with ID `61d5f311e722267a17c521f6`.
  479.40 KiB / 479.40 KiB ┃▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓┃100.00% 1.70 MiB/s 0s
'success'   '61d5f311e722267a17c521f6'      'G20502.22Rv1.2_fastqc.b64html` '490910'

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```








