

document.querySelector('title').textContent = 'Datalad on Helix';
Datalad on Helix


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Example](#int)
 |



Datalad provides users with access to a data sharing platform. Users can upload, version and share datasets as well as download datasets using datalad.



### Web site


* [Home page](https://www.datalad.org/)
* [Github page](https://github.com/datalad)


Documentation
* [Datalad Documentation](http://handbook.datalad.org/en/latest/)


Important Notes
* Please use datalad on an **Helix** and not on Biowulf.



Example
Establish a connection to Helix and run the program to download selected datasets from neurovault.   
Sample session (user input in **bold**):



```

[user@helix]$ **module load datalad git**
[+] Loading datalad  0.13.0rc2  on helix 
[+] Loading singularity  3.6.4  on helix
[+] Loading git 2.29.2  ... 

[user@helix]$ **cd /data/$USER** 

[user@helix]$ **datalad install -r ///neurovault**
install(ok): /data/user/neurovault (dataset)                                                                            
[INFO   ] Installing Dataset(/data/user/neurovault) to get /data/user/neurovault recursively            
install(ok): /data/user/neurovault/snapshots (dataset)                                                                  
action summary:                                                                                                                         
  install (ok: 2)

[user@helix]$ **cd neurovault**

[user@helix]$ **datalad get snapshots/1003/**
get(ok): snapshots/1003/13873.nii.gz (file) [from origin...]                                                                            
get(ok): snapshots/1003/13874.nii.gz (file) [from origin...]
get(ok): snapshots/1003/13875.nii.gz (file) [from origin...]
get(ok): snapshots/1003/13876.nii.gz (file) [from origin...]
get(ok): snapshots/1003/13877.nii.gz (file) [from origin...]
get(ok): snapshots/1003/13878.nii.gz (file) [from origin...]
get(ok): snapshots/1003/13879.nii.gz (file) [from origin...]
get(ok): snapshots/1003/13880.nii.gz (file) [from origin...]
get(ok): snapshots/1003/13881.nii.gz (file) [from origin...]
get(ok): snapshots/1003/13882.nii.gz (file) [from origin...]
  [10 similar messages have been suppressed]
get(ok): snapshots/1003 (directory)
action summary:
  get (notneeded: 1, ok: 21)


```





