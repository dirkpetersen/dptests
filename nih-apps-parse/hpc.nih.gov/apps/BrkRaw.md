

document.querySelector('title').textContent = 'BrkRaw: A comprehensive tool to access raw Bruker Biospin MRI data.';
BrkRaw: A comprehensive tool to access raw Bruker Biospin MRI data.


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
 |



The ‘BrkRaw’ is a python module designed to provide a comprehensive tool to access raw data acquired from Bruker Biospin preclinical MRI scanner. This module is also compatible with the zip compressed data to enable use of the archived data directly.
The module is comprised of four components, including graphical user interface (GUI), command-line tools, high-level and low-level python APIs.



### References:


* Lee, Sung-Ho, Ban, Woomi, & Shih, Yen-Yu Ian.   

 *BrkRaw/bruker: BrkRaw v0.3.3 (Version 0.3.3). Zenodo. June 4, 2020.*  

[DOI:](http://doi.org/10.5281/zenodo.3877179)  http://doi.org/10.5281/zenodo.3877179


Documentation
* [BrkRaw Github page](https://github.com/BrkRaw/bruker)


Important Notes
* Module Name: BrkRaw (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **BRKRAW\_HOME**  installation directory
	+ **BRKRAW\_BIN**       executable directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive** 
[user@cn4471 ~]$ **module load BrkRaw** 
[+] Loading BrkRaw 0.3.6  ...

```

Copy sample data to the current folder:

```

[user@cn4471 ~]$ **brkraw -h**
usage: brkraw [-h] [-v] command ...

BrkRaw command-line interface

optional arguments:
  -h, --help     show this help message and exit
  -v, --version  show program's version number and exit

Sub-commands:
  To run this command, you must specify one of the functions listedbelow
  next to the command. For more information on each function, use -h next to
  the function name to call help document.

  command        description
    info         Prints out the information of the internal contents in Bruker
                 raw data
    gui          Run GUI mode
    tonii        Convert a single raw Bruker data into NifTi file(s)
    tonii_all    Convert All raw Bruker data located in the input directory
    bids_helper  Creates a BIDS datasheet for guiding BIDS data converting.
    bids_convert
                 Convert ALL raw Bruker data located in the input directory
                 based on the BIDS datasheet

```

End the interactive session:

```

[user@cn4471 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





