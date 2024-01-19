

document.querySelector('title').textContent = 'Parmed on Biowulf';
Parmed on Biowulf


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



ParmEd is a package designed to facilitate creating and easily manipulating molecular systems that are fully described by a common classical force field. Supported force fields include Amber, CHARMM, AMOEBA, and several others that share a similar functional form (e.g., GROMOS).

ParmEd is capable of reading and writing to a wide array of different file formats, like the Amber topology and coordinate files, CHARMM PSF, parameter, topology, and coordinate files, Tinker parameter, topology, and coordinate files, and many others. The expressive central data structure (the Structure class) makes it easy to quickly and safely manipulate a chemical system, its underlying topology, and force field parameters describing its potential energy function.

There are two parts of ParmEd---a documented API that you can incorporate into your own Python scripts and programs, and a GUI/CLI pair of programs that provide a means to quickly perform various modifications to chemical systems for rapid prototyping.

The API also provides bindings to the OpenMM library, permitting you to carry out full molecular dynamics investigations using ParmEd on high-performant hardware, like AMD and NVidia GPUs.



Documentation
* [ParmEd site](https://github.com/ParmEd/ParmEd) on github


Important Notes
This application requires a [graphical connection using NX](/docs/connect.html#nx)


* Module Name: parmed (see [the modules page](/apps/modules.html) for more information)



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load parmed**

[user@cn3144 ~]$  **parmed --help**
parmed --help
usage: parmed [-h] [-v] [-i FILE] [-p <prmtop>] [-c <inpcrd>] [-O] [-l FILE] [--prompt PROMPT] [-n] [-e] [-s] [-r]
              [<prmtop>] [<script>]

positional arguments:
  <prmtop>              Topology file to analyze.
  <script>              File with a series of ParmEd commands to execute.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

Input Files:
  -i FILE, --input FILE
                        Script with ParmEd commands to execute. Default reads from stdin. Can be specified multiple times to process
                        multiple input files.
  -p <prmtop>, --parm <prmtop>
                        List of topology files to load into ParmEd. Can be specified multiple times to process multiple topologies.
  -c <inpcrd>, --inpcrd <inpcrd>
                        List of inpcrd files to load into ParmEd. They are paired with the topology files in the same order that each
                        set of files is specified on the command-line.

Output Files:
  -O, --overwrite       Allow ParmEd to overwrite existing files.
  -l FILE, --logfile FILE
                        Log file with every command executed during an interactive ParmEd session. Default is parmed.log

Interpreter Options:
  These options affect how the ParmEd interpreter behaves in certain cases.

  --prompt PROMPT       String to use as a command prompt.
  -n, --no-splash       Prevent printing the greeting logo.
  -e, --enable-interpreter
                        Allow arbitrary single Python commands or blocks of Python code to be run. By default Python commands will not
                        be run as a safeguard for your system. Make sure you trust the source of the ParmEd command before turning this
                        option on.

Error Handling:
  These options control how ParmEd handles various errors and warnings that appear occur during the course of Action execution

  -s, --strict          Prevent scripts from running past unrecognized input and actions that end with an error. In interactive mode,
                        actions with unrecognized inputs and failed actions prevent any changes from being made to the topology, but
                        does not quit the interpreter. This is the default behavior.
  -r, --relaxed         Scripts ignore unrecognized input and simply skip over failed actions, executing the rest of the script.
                        Unrecognized input in the interactive interpreter emits a non-fatal warning.
[user@cn4272 ~]$

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. parmed.sh). For example:



```

#!/bin/bash
set -e
module load parmed
parmed .....

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] parmed.sh
```













