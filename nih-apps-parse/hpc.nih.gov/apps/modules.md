

document.querySelector('title').textContent = ' Environment Modules on Biowulf & Helix';

Environment Modules on Biowulf & Helix



|  |
| --- |
| 
Quick Links
[Video Tutorial](#vid)
[Available Modules](#avail)
[Load a Module](#load)
[Unload a Module](#unload)
[See Loaded Modules](#view)
[Switch Application Version](#switch)
[Examine a modulefile](#examine)
[Personal Modulefiles](#personal)
[Using modules in scripts](#scripts)
[Nested Modules](modules_advanced.html#nested)
[Shared Modules](#shared)
[Output from modules](#output)
 |



Helix and Biowulf now use environment modules to dynamically set up environments for different applications. Users can type commands like:

> 
> ```
> 
> module load ThisApp
> module unload ThisApp
> 
> ```
> 



to set up the environment for a particular application. The 'module load' command will set PATH, LD\_LIBRARY\_PATH and other environment 
variables as necessary. Any previous environment setups in your .bashrc or .cshrc files will continue to work, but modules are generally easier to use.

The Biowulf application pages have been updated to describe the modules command for each application. On 4 June 2018, Helix transition to an ['interactive data transfer node'](https://hpc.nih.gov/docs/rhel7.html#helix) and a limited number of apps related to data transfer are available on Helix. This page contains some general information about modules.

We use the [Lmod Environment Modules system](https://www.tacc.utexas.edu/tacc-projects/lmod), developed at the TACC. 



  


Video Tutorial



What modules are available?

The command module avail will show all available modules.


```

[biowulf ~]$ **module avail**
------------------------- /usr/local/lmod/modulefiles --------------------------
   GSL/1.16                    mpich2/1.5/gnu
   ImageMagick/6.8.9           oncotator/-1.2.8.1
   cmake/3.0.0                 oncotator/0.1.2.8.0
   dot                         oncotator/1.2.8.0
   fcgene/1.0.6                oncotator/1.2.8.1       (D)
   fcgene/1.0.7         (D)    openmpi/1.8.1/gnu-eth   (D)
   fftw/3.3.4/gnu       (D)    openmpi/1.8.1/intel-eth
   fftw/3.3.4/intel            openmpi/1.8.1/pgi-eth
   fftw/3.3.4/pgi              parallel/20140722
   graphviz/2.34               pgi/14.7
   hdf5/1.8.13                 plink/1.06
   inspect/2012.01.09          plink/1.07              (D)
   intel/2013_sp1.3.174        python/2.7.8
   java/1.7.0_25               testnest/1.0
   java/1.8.0_11        (D)    testnest/1.1            (D)
   mathematica/9.0             trinity/r20140413
   meme/4.9.1                  trinity/r20140717       (D)
   meme/4.10.0          (D)    use.own
   module-info                 vmd/1.9
   modules                     vmd/1.9.1               (D)
[...etc...]   

```


Since there are so many applications and versions installed on Biowulf, this list can be overwhelming. 
If you add the '-d' flag, only the default version of each application will be listed.

```

[biowulf ~]$ **module -d avail**
------------------------- /usr/local/lmod/modulefiles --------------------------
   GSL/1.16              intel/2013_sp1.3.174     parallel/20140722
   ImageMagick/6.8.9     java/1.8.0_11            pgi/14.7
   cmake/3.0.0           mathematica/9.0          plink/1.07
   dot                   meme/4.10.0              python/2.7.8
   fcgene/1.0.7          module-info              testnest/1.1
   fftw/3.3.4/gnu        modules                  trinity/r20140717
   graphviz/2.34         mpich2/1.5/gnu           use.own
   hdf5/1.8.13           oncotator/1.2.8.1        vmd/1.9.1
   inspect/2012.01.09    openmpi/1.8.1/gnu-eth
  [...etc...] 
   
```


If you prefer the output of 'module -d avail', you can set up an alias in your .bashrc file:

```

alias ml="module -d avail"

```


To see the versions available for a particular application, use module avail appname. e.g.


```

[biowulf ~]$ **module avail gromacs**

---------- /usr/local/lmod/modulefiles -----------
gromacs/4.0.4          gromacs/4.5.1          gromacs/4.5.5(default)
gromacs/4.0.7          gromacs/4.5.3          gromacs/4.6-dev

```

This command, along with other search commands (spider, list) supports the use of
regular expressions. This can be helful for packages with difficult to search names
such as 'R':



```

[biowulf ~] **module -r avail '^R$'**

------------------------- /usr/local/lmod/modulefiles --------------------------
   R/Rdev_gcc-4.9.1     R/3.3.0_gcc-4.7.4    R/3.4.0_gcc-4.9.1     (D)
   R/Rdev_gcc-6.2.0     R/3.3.0_gcc-4.9.1    R/3.4.0_gcc-6.2.0
   R/2.14_gcc-4.4.7     R/3.3.1_gcc-4.9.1    R/3.4.3_gcc-4.9.1
   R/3.2.3_gcc-4.9.1    R/3.3.1_gcc-6.2.0    R/3.5.0beta_gcc-4.9.1
   R/3.2.5_gcc-4.9.1    R/3.3.2_gcc-4.9.1

  Where:
   (D):  Default Module

Use "module spider" to find all possible modules.
Use "module keyword key1 key2 ..." to search for all possible modules matching
any of the "keys".

```


A very useful command is 'module spider', which will run a case-insensitive search for any module that includes that word. e.g.

```

[biowulf ~] **module spider bed**
--------------------------------------------------
  bedops: bedops/2.4.3
-------------------------------------------------
    This module can be loaded directly: module load bedops/2.4.3


-------------------------------------------------
  bedtools:
-------------------------------------------------
     Versions:
        bedtools/2.21.1
        bedtools/2.22.0

-------------------------------------------------
  To find detailed information about bedtools please enter the full name.
  For example:

     $ module spider bedtools/2.22.0
------------------------------------------------

```


Load a Module

To load a module, use the command module load appname. The default version will get loaded. If you want a particular version, use module load appname/version e.g


```

[biowulf ~]$ **module list**
No Modulefiles Currently Loaded.

[biowulf ~]$ **module load gromacs**

[biowulf ~]$ **module list**
Currently Loaded Modules:
  1) gromacs/4.5.5

```


```

[biowulf ~]$ **module load gromacs/4.5.3**

[biowulf ~]$ **module list**
Currently Loaded Modules:
  1) gromacs/4.5.3   

```


Multiple modules can be loaded in a single command. e.g.

```

[biowulf ~]$ **module load plinkseq macs bowtie**

[biowulf ~]$ **module list**
Currently Loaded Modules:
  1) plinkseq/0.08    2) macs/2.0.10      3) bowtie/2-2.0.2

```


Unload a Module

To unload a module, use the command module unload Appname. e.g.

```

[biowulf ~]$ **module list**
Currently Loaded Modules:
  1) gromacs/4.5.3
[biowulf ~]$ **module unload gromacs/4.5.3**

```


See what modules have been loaded

The command module list will show what modules have been loaded. 

```

[biowulf ~]$ **module list**
Currently Loaded Modules:
  1) qiime        2) perl/5.8.9

```


Switch to a different version of an application

The command module switch appname appname/version will switch from one version of an application to another. .


```

[biowulf ~]$ **module list**
Currently Loaded Modules:
  1) perl/5.8.9
[biowulf ~]$ **module switch perl perl/5.12.1**
[biowulf ~]$ **module list**
Currently Loaded Modules:
  1) perl/5.12.1

```


If you load a module, then load another version of the same module, the first one will be unloaded. For example:

```

[biowulf ~]$ **module load fcgene**
[biowulf ~]$ **module list**

Currently Loaded Modules:
  1) fcgene/1.0.7

[biowulf ~]$ **module load fcgene/1.0.6**

The following have been reloaded with a version change:
  1) fcgene/1.0.7 => fcgene/1.0.6

[biowulf ~]$ **module list**

Currently Loaded Modules:
  1) fcgene/1.0.6
  
```



Examine a modulefile

If you want to see what changes a module will make to your environment without loading it, use the 'module display' command. e.g.


```

[user@biowulf] **module display meme**

----------------------------------------------------------------------------------------------
   /usr/local/lmod/modulefiles/meme/4.10.0:
---------------------------------------------------------------------------------------------
whatis("Sets up MEME 4.10.0 - MPI version ")
prepend_path("PATH","/usr/local/OpenMPI/1.8.1/gnu/eth/bin")
prepend_path("LD_LIBRARY_PATH","/usr/local/OpenMPI/1.8.1/gnu/eth/lib")
prepend_path("PATH","/usr/local/apps/meme/4.10.0/bin")

```


The output is in the **lua** programming language, but is usually easily understandable. For a table of lua-tcl-bash equivalents, see 
[this page](modules_advanced.html).


If you prefer to set paths directly instead of using modules, you can use the information from "**module display**" to do so. 


Set up personal modulefiles

Create a directory called "**modulefiles**" in your home directory. Then type the command: "**module use --append ~/modulefiles**" and any personal module files in this directory will become available for you to list, load or unload. You can make copies of the system modulefiles in "**/usr/local/lmod/modulefiles**" to use as templates. 

```

[biowulf ~]$ **module use --append ~/modulefiles**

[biowulf ~]$ **module avail**

---------- /usr/local/lmod/modulefiles -----------
afni/18Oct2012(default) cuda/3.1                modules                 namd/2.9-ipath
afni/18Oct2012-openmp   cuda/4.0.17             mpich2-x86_64           null
afni/21Dec2011          dot                     namd/2.7+plumed1.3      perl/5.12.1
[...etc...]

---------- /home/user/modulefiles ----------
hmmer/2.3



```


As you see above, after adding '~/modulefiles' to the module path, your private modules are now displayed. If a personal module has the same name as a system module, the system module will have priority. You can also use **module use --prepend ~/modulefiles** to give your personal modules priority over system modules. 

Using modules in scripts

The '**module**' command can be used in any batch script or other script. A typical script on the Biowulf login node would be 


```

#!/bin/bash

module load somemodule
somecommands

```


Please note that, if your login shell is bash and you would like to submit a job using a script written in tcsh, incompatibilities between bash and tcsh will prevent the export of parts of your environment, such as defined functions. This means that the module command will not work automatically. To get around this, you should load the module setup script for csh as part of your submission script:


```

#!/bin/tcsh

source /etc/profile.d/modules.csh
module load somemodule
somecommands

```


A typical Biowulf batch script to be submitted to the Biowulf batch system would be:

```

#!/bin/bash

module load afni
[...some afni command...]

```


Likewise, Biowulf swarm jobs could either include the module command in each line, or be submitted with --module somemodule. See the [Biowulf swarm page](/apps/swarm.html) for more information.

[Working with Modules within Perl and Python](http://www.nersc.gov/users/computational-systems/genepool/user-environment/working-with-modules-within-perl-and-python/), at the NERSC site.


For information about Nested Modules (a modulefile that loads and unloads other modules), see 
[this page](modules_advanced.html)


Shared Modules

To set up modules that will be shared by a group, the easiest way is to create a 'modulefiles' directory in some area that is shared by the group, and put the group modules into that directory. Each user in the group who wants to access these modules can then create a link from the shared area to their own /home/user/modulefiles/ area. 
For example, suppose the group has a shared area called /data/DBCImaging.

```

biowulf% **mkdir /data/DBCImaging/modulefiles**

```

Create the desired modulefiles in this directory. If a personal/group 
module has the same name as a system module, the system module will get lower priority, and your personal module will be loaded first.

Now each user in the group who wants to use these shared modules should create a link from their /home/user/modulefiles area.

```

biowulf% **ln -s /data/DCBImaging/modulefiles /home/user/modulefiles/shared**
biowulf% **module load use.own**
biowulf% **module avail**

----------------------- /home/user/modulefiles -------------------------------
blast                       shared/bowtie/0.12.8        shared/bowtie/2-2.0.6
freesurfer/5.1              shared/bowtie/0.12.9        shared/bowtie/2-2.1.0
myngs                       shared/bowtie/2-2.0.0-beta7 shared/octave_test/3.2.4
nested                      shared/bowtie/2-2.0.2       shared/octave_test/3.6.1
shared/bowtie/0.12.5        shared/bowtie/2-2.0.5       susantest

---------------------------- /usr/local/lmod/modulefiles ---------------
BEAST/1.5.2                    eman2/2.06                     openmpi/1.6/gnu/eth
BEAST/1.5.4                    eman2/2.06_121120              openmpi/1.6/gnu/ib
BEAST/1.6.1                    eman2/2.06_ib                  openmpi/1.6/gnu/ipath
BEAST/1.6.2                    express/1.2.2                  openmpi/1.6/intel/eth
BEAST/1.7.1                    express/1.3.1(default)         openmpi/1.6/intel/ib
[...]

biowulf% **module load shared/octave\_test/3.2.4**

biowulf% **module list**
Currently Loaded Modules:
  1) use.own                    2) shared/octave_test/3.6.1

```


Users who do not have access to the shared area (in this case /data/DBCImaging) will not be able to see the shared modules.


Output from Modules

Most modules will silently set up paths and environments. Some modules need to set up additional dependencies apart from the main 
program (e.g. a scientific program that requires a particular version of Python will load that Python version into your path as well), and 
in that case, we have set up the module to print the additional dependencies that are being loaded. For example:


```

[biowulf ~]$ module load sratoolkit
java/1.7.0 is loaded

```


These messages are useful for interactive use, but may be annoying or problematic when running large-scale batch jobs. In that case, to 
get rid of the messages, you should use 'module -q load' (-q = quiet)


```

#!/bin/bash

module -q load somemodule

program [options] [arguments]....

```


Documentation

[Lmod: A new environment module system](https://lmod.readthedocs.io/)
































































































