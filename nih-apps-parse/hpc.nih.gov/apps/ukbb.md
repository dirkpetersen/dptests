

document.querySelector('title').textContent = 'UKBB on Helix';
UKBB on Helix


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Basic Usage](#int)
 |



UKBB are the utilities to download data from the UK Biobank project.



### Web site


* [Home page](https://biobank.ctsu.ox.ac.uk/)


Documentation
* [UKBB Documentation](http://biobank.ctsu.ox.ac.uk/showcase/)


Important Notes
* Module Name: ukbb (see [the modules page](/apps/modules.html) for more information)
* UKBB utitlies are available only on Helix.



Basic usage
  
Sample session (user input in **bold**):



```

[user@helix ~]$ **module load ukbb**

[user@helix ~]$ **ukbfetch**

ukbfetch on unx - ver Mar 14 2018 14:21:29 - using Glibc2.12(stable)
Run start : 2018-10-12T14:05:06 
Must specify encoded_id for participant (-e flag)
Usage: ukbfetch parameters...
 -a	authentication file (application_id + 24-char key)
 -b	batch file containing list of participants and datafiles
 -d	datafile name
 -e	encoded id for participant
 -h	show this usage message then exit
 -i	show program version information only then exit
 -m	maximum datafiles to fetch (batch mode only, capped at 500)
 -o	name of output file recording successful fetches
 -s	starting line (batch mode only)
 -v	verbose mode on

Compiled : Mar 14 2018 14:21:29

[user@helix ~]$

```








