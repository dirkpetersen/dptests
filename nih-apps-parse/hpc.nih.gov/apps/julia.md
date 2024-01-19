

document.querySelector('title').textContent = 'JULIA on Biowulf';
JULIA on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Common pitfalls](#gotcha)
[Interactive job](#int) 
[Package management](#pkg) 
[Julia + Jupyter notebooks](#jupyter)
[Pluto notebooks](#pluto)
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
 |



[![julia logo](/images/julia.png)](http://julialang.org)

Julia is a flexible dynamic language for scientific and numerical computing. Julia's performance is comparable to traditional statically-typed languages.




Julia has a high level syntax and has good support for interactive use. Julia is provided under the MIT license (free) and the source code is publically available on Github. 



 

Documentation
* [Julia Main Site](https://julialang.org)
* [Julia Documentation](https://docs.julialang.org)
* [Julia GitHub page](https://github.com/JuliaLang)



Notes
* Module Name: julialang (see [the modules page](/apps/modules.html) for more information)




Common pitfalls


Multithreading
Some Julia libraries will attempt to use all CPUs on a compute node. Unless all CPUs have been allocated
this will result in overloading the job. Please be sure to set `OMP_NUM_THREADS` 
and `OPENBLAS_NUM_THREADS` explicitly for all code that can potentially multithread. If using the 
`Distributed`  package please first set `OPENBLAS_NUM_THREADS` and `OMP_NUM_THREADS` to 1
to prevent multi-threading, and then set the number of cores in your Julia script as follows: 


```

# use distributed package to run parallel julia scripts
using Distributed
# launch workers to match the allocated CPUs
num_cores = parse(Int, ENV["SLURM_CPUS_PER_TASK"])
addprocs(num_cores)

```




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

```

There may be multiple versions of julia available. An easy way of selecting the
version is to use [modules](/apps/modules.html).To see the modules
available, type



```

[user@cn3144 ~]$ **module avail julia**
--------------------------------------- /usr/local/lmod/modulefiles ----------------------------------------
   julia/0.6.2    julia/1.3.1    julia/1.6.0        julialang/0.6.2    julialang/1.3.1    julialang/1.6.0
   julia/1.0.0    julia/1.4.0    julia/1.6.1 (D)    julialang/1.0.0    julialang/1.4.0    julialang/1.6.1 (D)
   julia/1.1.0    julia/1.5.0    julia/1.7.1        julialang/1.1.0    julialang/1.5.0    julialang/1.7.1


  Where:
   L:  Module is loaded
   D:  Default Module

Use "module spider" to find all possible modules.
Use "module keyword key1 key2 ..." to search for all possible modules matching any of the "keys".

```

Note that the julia and julialang modules are identical to accomodate both conventional
names of the language. You can now load the module that corresponds to the desired Julia version and start an 
interactive Julia session as follows:



```

[user@cn3144 ~]$ **module load julialang/1.7.1**
[+] Loading git 2.34.1  ... 
[+] Loading julialang  1.7.1

[user@cn3144 ~]$ **julia**
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.7.1 (2021-12-22)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia>

julia> **5 / 2**
2.5

julia> **exit()**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Package management
To install, update, and remove packages you can use Julia's package manager, [Pkg](https://docs.julialang.org/en/latest/stdlib/Pkg/). Pkg comes with a REPL shell that you can start from within Julia by typing **]**. To exit the REPL enter **^C** or backspace. Within the REPL you can obtain a list of available commands by typing **?**:



```

julia> **]**
(v1.7) pkg> **?**
  Welcome to the Pkg REPL-mode. To return to the julia> prompt, either press backspace when the input line is empty or press Ctrl+C.
  Full documentation available at https://pkgdocs.julialang.org/

  Synopsis

  pkg> cmd [opts] [args]

  Multiple commands can be given on the same line by interleaving a ; between the commands. Some commands have an alias, indicated below.

  Commands

  activate: set the primary environment the package manager manipulates
  add: add packages to project
  build: run the build script for packages
  develop, dev: clone the full package repo locally for development
  free: undoes a pin, develop, or stops tracking a repo
  gc: garbage collect packages not used for a significant time
  generate: generate files for a new project
  help, ?: show this message
  instantiate: downloads all the dependencies for the project
  pin: pins the version of packages
  precompile: precompile all the project dependencies
  redo: redo the latest change to the active project
  remove, rm: remove packages from project or manifest
  resolve: resolves to update the manifest from changes in dependencies of developed packages
  status, st: summarize contents of and changes to environment
  test: run tests for packages
  undo: undo the latest change to the active project
  update, up: update packages in manifest
  registry add: add package registries
  registry remove, rm: remove package registries
  registry status, st: information about installed registries
  registry update, up: update package registries


```

To obtain a list of installed packages, for example, type **status**:



```

(v1.7) pkg> **status**
    Status `.julia/environments/v1.7/Project.toml` (empty project)

```

To install a package use **add**, for example:



```

(v1.7) pkg> **add Flux**
  Updating registry at `~/.julia/registries/General`
  Updating git-repo `https://github.com/JuliaRegistries/General.git`
   Resolving package versions...
   Installed IrrationalConstants ─ v0.1.1
   Installed RealDot ───────────── v0.1.0
   Installed Adapt ─────────────── v3.3.3
   Installed DiffRules ─────────── v1.9.0
   [...]
   [8e850b90] + libblastrampoline_jll
   [8e850ede] + nghttp2_jll
   [3f19e933] + p7zip_jll
Precompiling project...
   69 dependencies successfully precompiled in 95 seconds

```

To update your packages use **up**:



```

(v1.7) pkg> **up**
    Updating registry at `~/.julia/registries/General`
    Updating git-repo `https://github.com/JuliaRegistries/General.git`
  No Changes to `/spin1/home/linux/$USER/.julia/environments/v1.7/Project.toml`
  No Changes to `/spin1/home/linux/$USER/.julia/environments/v1.7/Manifest.toml`
[ Info: We haven't cleaned this depot up for a bit, running Pkg.gc()...
      Active manifest files: 8 found
      Active artifact files: 15 found
      Active scratchspaces: 5 found
     Deleted no artifacts, repos, packages or scratchspaces

```

You can also use Pkg's REPL to create new environments and organize dependencies for each environment. For example the following creates a new environment, installs the Plots dependency in that environment, then switches back to the v1.1 environment:



```

(v1.7) pkg> **activate ~/.julia/environments/new\_environment/**
[ Info: activating new environment at ~/.julia/environments/new_environment.

(new_environment) pkg> **status**
    Status `~/.julia/environments/new_environment/Project.toml` (empty environment)

(new_environment) pkg> **add Plots**
   Installed Xorg_xcb_util_renderutil_jll ─ v0.3.9+1
   Installed JpegTurbo_jll ──────────────── v2.1.0+0
   Installed FFMPEG ─────────────────────── v0.4.1
   [...]
  [8e850ede] + nghttp2_jll
  [3f19e933] + p7zip_jll
    Building GR → `~/.julia/scratchspaces/44cfe95a-1eb2-52ea-b672-e2afdf69b78f/4a740db447aae0fbeb3ee730de1afbb14ac798a1/build.log`
Precompiling project...


(new_environment) pkg> **status**
    Status `~/.julia/environments/new_environment/Project.toml`
  [91a5bcdd] Plots v1.25.6

(new_environment) pkg> **activate ~/.julia/environments/v1.7/**

(v1.7) pkg> 


```

For the full documentation of Pkg's REPL, check the [Pkg docs](https://docs.julialang.org/en/latest/stdlib/Pkg/). There is a searchable list of registered Julia packages at [Julia packages](https://pkg.julialang.org/). 

There are two ways to use this project/environment in a batch script non interactively:



```

#! /bin/bash

module load julia
julia --project=$HOME/.julia/environments/new_environment my_simulation_of_everything.jl

```

or



```

#! /bin/bash

module load julia
export JULIA_PROJECT==$HOME/.julia/environments/new_environment
julia my_simulation_of_everything.jl

```

### Reducing the size of the ~/.julia folder


The size of the `$HOME/.julia` folder where packages, environments, etc
are stored can grow quite large. Here are some ways to reduce size:



Removing environments
There isn't an obvious way to delete an environment from the
 julia Pkg mode. To delete an environment and all its packages
 first delete the environment folder from the file system and then run
 a Pkg garbage collection to remove the actual packages. Using the example above:
 
```

$ **rm -rf ~/.julia/new\_environment**
$ **julia**
julia> **]**
(v1.7) pkg> **gc --all**
     
```


Registries
julia 1.7 changed the way registries are handled. Once you update to 1.7
 it may be a good idea to remove and then add back registries.
 
```

(v1.7) pkg> **registry status**
Registry Status
 [23338594] General
(v1.7) pkg> **registry rm General**
    Removing registry `General` from ~/.julia/registries/General
(v1.7) pkg> **registry add General**
    
```

Relocate .julia folder
If the folder is still too large it can be relocated to /data.
 
```

$ **cd ~**
$ **mv .julia /data/$USER/.julia**
$ **ln -s /data/$USER/.julia**
    
```




Julia on Jupyter
You can have access to a Julia kernel through Jupyter by installing the IJulia package in your home directory and setting up the 
appropriate tunneling (as described in [Jupyter on Biowulf](jupyter.html)). 
To install the IJulia package do the following on a Julia interactive session first start the pkg mode with **]**



```

(v1.7) pkg> **add IJulia**
  Updating registry at `~/.julia/registries/General`
  Updating git-repo `https://github.com/JuliaRegistries/General.git`
 Resolving package versions...
 Installed Rmath_jll  v0.2.2+0
 Installed Gettext_jll  v0.20.1+1
 [...]
(v1.7) pkg>


```

In addition to installing the IJulia package, this will also install a jupyter kernel
in your home directory under `~/.local/share/jupyter/kernels` which will then
be included in the kernels menu of jupyter notbook and jupyter lab.


Note that we do not install Julia packages centrally so this has to be done by each
user separately



Pluto notebooks
Pluto notebooks can be setup similarly to Jupyter notebooks. Allocate an interactive session,
start julia, enter the pkg mode, and install the Pluto package. Note that the `--tunnel`
sets up a tunnel between the login node and the compute node. See also our [tunneling docs](/docs/tunneling/).



```

[user@biowulf]$ **sinteractive --cpus-per-task=4 --tunnel**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job
[user@cn3144 ~]$ **module load julialang/1.7.1**
[+] Loading git 2.34.1  ... 
[+] Loading julialang  1.7.1

[user@cn3144 ~]$ **julia**
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.7.1 (2021-12-22)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia>
julia>> **]**
(v1.7) pkg> **add Pluto**
    Updating registry at `~/.julia/registries/General.toml`
   Resolving package versions...
   Installed IniFile ──────────── v0.5.1
   Installed RelocatableFolders ─ v0.2.0
   [...]
(v1.7) pkg> ^C
julia>

```

Then start the notebook



```

julia> **using Pluto**
julia> # the following should not cause any errors unless you forgot to use --tunnel
          # when starting the sinteractive session
julia> **port = parse(Int, get(ENV, "PORT1", ""))**
julia> **Pluto.run(host="127.0.0.1", port=port, launch\_browser=false)**

Go to http://localhost:34090/?secret=XXXXXX in your browser to start writing ~ have fun!

Press Ctrl+C in this terminal to stop Pluto

    Updating registry at `~/.julia/registries/General.toml`
    Updating registry done ✓


```

At this point you should set up a tunnel from your local computer to the login
node according to the instructions printed out by sinteractive (
see also [Jupyter on Biowulf](jupyter.html)). One that is
set up you can visit the URL printed by Pluto from your local broswer.

Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. julia.sh). For example:



```

#!/bin/bash
set -e
module load julialang
myscript.jl < data.in > result.out

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] julia.sh
```


Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. julia.swarm). For example:



```

myscript.jl < dataset1.in > result1.out
myscript.jl < dataset2.in > result2.out
myscript.jl < dataset3.in > result3.out
myscript.jl < dataset4.in > result4.out

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f julia.swarm [-g #] [-t #] --module julialang
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module julialang Loads the julialang module for each subjob in the swarm 
 | |
 | |
 | |














