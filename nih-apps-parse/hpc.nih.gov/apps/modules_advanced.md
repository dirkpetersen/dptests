

document.querySelector('title').textContent = ' Environment Modules on Biowulf & Helix Advanced Topics';

Environment Modules on Biowulf & Helix  

Advanced Topics

Examine a modulefile, contd.

The Module Environment System, "**Lmod**", currently in use on **biowulf** and **helix**, always reports the output of the "**module display**" command in the "**lua**" programming language, as shown in the above example. If you are more familiar with the "**tcl**" programming language, use **Table 1.** below to translate the "**Lmod**" output.



Table 1. Translations of some modulefile commands
| tcl command | lua command | shell | shell equivalent |
| prepend\_path PATH NewPath | prepend\_path("PATH","NewPath") | any | PATH=NewPath:$PATH |
| append\_path PATH NewPath | append\_path("PATH","NewPath") | any | PATH=$PATH:NewPath |
| setenv VarName varvalue | setenv(VarName,varvalue) | bashksh | export VarName=varvalue |
| csh | setenv VarName varvalue |
| module-whatis "Whatis message" | whatis("Whatis message") |  |



If you prefer to set paths directly instead of using modules, you can use the information from "**module display**" to do so. 


Nested modules

You might want to set up a personal modulefile that loads and unloads a collection of other modules. The advantage to writing a modulefile in this way is: if a module "**A**" requires modules "**B**","**C**", and "**D**" (or more!), the single module command,
"**module load A**" will load module "**A**" and all its required modules, thus avoiding the commands
"**module load B"**, "**module load C**", "**module load D**" (or more!).

Below is an example of such a module, called '**myngs**', which loads a specific version of **Novocraft**, **Cufflinks**, and **Tophat**. You will need to put this module file into "**/home/user/privatemodules/**".


```

#%Module1.0
 
# this module is called myngs
  
module-whatis "Loads  novocraft/3.00.02 cufflinks/2.0.1   tophat/2.0.8/bowtie2-2.1.0"
 
set mode [module-info mode]
#
if { $mode == "whatis" } {  exit }
#
foreach modulename { novocraft/3.00.02 cufflinks/2.0.1   tophat/2.0.8/bowtie2-2.1.0 } {
  module load $modulename
  if { $mode == "load" }   { puts stderr "$modulename is loaded" }
  if { $mode == "unload" } { puts stderr "$modulename is unloaded" }
}

```


The example module "**myngs**" has two notable features:
1. For each required module, with name "**$modulename**" used as an item in the **foreach** loop,   

there is a command "**module load $modulename**" but no "**module unload $modulename**". 

This works only because the new Environment Module manager, "**lmod**", has this behavior:

	* When the user performs "**module load myngs**" in a shell, each "**module load $modulename**" in the **modulefile** performs a "**load**" operation.
	* When the user performs "**module unload myngs**" in a shell, each "**module load $modulename**" in the **modulefile** performs an "**unload**" operation.

- The lines with "**puts stderr**" make a report only as a reminder, and are optional.



Here is a sample session using the above module.

```

[user@helix ~]$ **module load use.own**

[user@helix ~]$ **module list**
Currently Loaded Modules:
  1) use.own

[user@helix ~]$ **module load myngs**
novocraft/3.00.02 is loaded
cufflinks/2.0.1 is loaded
tophat/2.0.8/bowtie2-2.1.0 is loaded

[user@helix ~]$ **module list**
Currently Loaded Modules:
  1) use.own                      3) cufflinks/2.0.1              5) myngs
  2) novocraft/3.00.02            4) tophat/2.0.8/bowtie2-2.1.0

[user@helix ~]$ **module unload myngs**
novocraft/3.00.02 is unloaded
cufflinks/2.0.1 is unloaded
tophat/2.0.8/bowtie2-2.1.0 is unloaded

[user@helix ~]$ **module list**
Currently Loaded Modules:
  1) use.own

```


Tip on Best-Practices When Writing Personal Modulefiles

One purpose of a **modulefile** for an application is to adjust environment variables to the needs of the application.
Particularly important are variables which have a value specifying  
 a "search path" — a colon-separated list of directory pathnames in search order.

 For example:



|  |  |
| --- | --- |
| Variable Name | Description of Directory Contents |
| **PATH** | Executable programs |
| **MANPATH** | Linux man pages |
| **LD\_LIBRARY\_PATH** | Shared libraries |
| **PERL5LIB** | Perl modules |
| **PYTHONPATH** | Python packages |



The safest way to set these environment variables is to follow the rules:
* A **modulefile** performs a "**module load** *\*REQUIRED\**" operation  

 for each of its "*\*REQUIRED\**" modules

* Each **modulefile** adjusts the environment variables for its own application,  

and never those for its "*\*REQUIRED\**" modules.




































