

document.querySelector('title').textContent = 'SPHIRE: semi-automated processing of single particle electron cryo-microscopy (cryo-EM) data ';
**SPHIRE: semi-automated processing of single particle electron cryo-microscopy (cryo-EM) data** 


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



SPHIRE (SPARX for High-Resolution Electron Microscopy) is an open-source, 
user-friendly software suite for the semi-automated
processing of single particle electron cryo-microscopy (cryo-EM) data.
It allows fast and reproducible structure determination from cryo-EM images.



### References:


* T.Moriya, M.Saur, M.Stabrin, F.Merino, H.Voicu, Z.Huang, P.A.Penczek, S.Raunser, C.Gatsogiannis   

 *High-resolution Single Particle Analysis from Electron Cryo-microscopy
Images Using SPHIRE*  

 J. Vis Exp., 2017, **123**(21), p.55448; doi:10.3791/55448.


Documentation
* [SPHIRE Wiki page](http://sphire.mpg.de/wiki/doku.php)


Important Notes
* Module Name: SPHIRE (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **SPHIRE\_HOME**  installation directory
	+ **SPHIRE\_BIN**       executable directory
	+ **SPHIRE\_TEST**  test scripts



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=16g --gres=gpu:v100:1**
[user@@cn4471 ~]$ **module load SPHIRE**
[+] Loading SPHIRE  1.1

```

Run basic commands for the spire executable:

```

[user@cn4471 SPHIRE]**sphire -h**
Usage: sphire 
	The main SPHIRE GUI application. It is designed as the command generator for the SPHIRE single particle analysis pipeline.

Options:
  --version   show program's version number and exit
  -h, --help  show this help message and exit

[user@cn4471 SPHIRE]**sphire**
Python 2.7.14 | packaged by conda-forge | (default, Dec  9 2017, 16:18:43) 
Type "copyright", "credits" or "license" for more information.

IPython 5.8.0 -- An enhanced Interactive Python.
?         -> Introduction and overview of IPython's features.
%quickref -> Quick reference.
help      -> Python's own help system.
object?   -> Details about 'object', use 'object??' for extra details.
Welcome to the interactive SPARX-NoGUI Python interface, provided by ipython
   SPARX v4.0 (GITHUB: 2018-08-20 20:32)

 **In [1]: import mpi**

 **In [2]: Util.version()** 
   Source modification date: August 2018 release 

**In [3]: quit** 
[user@cn4471 SPHIRE]$ **e2version.py** 
EMAN 2.22 final (GITHUB: 2018-08-20 20:32 - commit: f4f3952 )
Your EMAN2 is running on: Linux-3.10.0-862.14.4.el7.x86_64-x86_64-with-centos-7.5.1804-Core 3.10.0-862.14.4.el7.x86_64
Your Python version is: 2.7.14

```

Copy test scripts to your current folder:

```

[user@@cn4471 ~]$ **cp -r $SPHIRE\_TEST/\* .**

```

Run a test:

```

[user@@cn4471 ~]$ **./test\_cmp.py** 
test_DotCmp (__main__.TestCmp)
test DotCmp ...................................... ... ok
test_FRCCmp (__main__.TestCmp)
test FRCCmp ...................................... ... ok
test_OptVarianceCmp (__main__.TestCmp)
test OptVarianceCmp .............................. ... ok
test_PhaseCmp (__main__.TestCmp)
test PhaseCmp .................................... ... ok
test_QuadMinDotCmp (__main__.TestCmp)
test QuadMinDotCmp ............................... ... ok
test_SqEuclideanCmp (__main__.TestCmp)
test SqEuclideanCmp .............................. ... ok
test_basic_cmp (__main__.TestCmp)
test basic cmp ................................... ... ok
test_variance (__main__.TestCmp)
test variance .................................... ... ok

----------------------------------------------------------------------
Ran 8 tests in 0.041s

OK
[user@cn4471 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





