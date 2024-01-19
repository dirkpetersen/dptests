

document.querySelector('title').textContent = 'ROBOT: a tool for working with Open Biomedical Ontologies.';
**ROBOT: a tool for working with Open Biomedical Ontologies.** 


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



ROBOT is a command-line tool and library for automating ontology development tasks, 
with a focus on Open Biological and Biomedical Ontologies (OBO).
It can be used as a command-line tool or as a library for any language on the Java Virtual Machine.



### References:


* R.C. Jackson, J.P. Balhoff, E. Douglass, N.L. Harris, C.J. Mungall, and J.A. Overton   

*ROBOT: A tool for automating ontology workflows,*    

[BMC Bioinformatics, vol. 20, July 2019.](https://link.springer.com/article/10.1186/s12859-019-3002-3)


Documentation
* [ROBOT GitHub page](https://github.com/ontodev/robot)
* [ROBOT Home page](http://robot.obolibrary.org)


Important Notes
* Module Name: robot (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **ROBOT\_HOME**  installation directory
	+ **ROBOT\_BIN**       executable directory
	+ **ROBOT\_SRC**       source code directory
	+ **ROBOT\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=4g**
[user@cn0911 ~]$**module load robot** 
[+] Loading singularity  3.10.5  on cn0833
[+] Loading robot  1.9.4
[user@cn0911 ~]$**robot**
usage: robot [command] [options] 
 --add-prefix  add prefix 'foo: http://bar' to the output
 --add-prefixes  add JSON-LD prefixes to the output
 --catalog  use catalog from provided file
 -h,--help print usage information
 --noprefixes do not use default prefixes
 -p,--prefix  add a prefix 'foo: http://bar'
 -P,--prefixes  use prefixes from JSON-LD file
 --strict use strict parsing when loading an ontology
 -V,--version print version information
 -v,--verbose increased logging
 -vv,--very-verbose high logging
 -vvv,--very-very-verbose maximum logging, including stack traces
 -x,--xml-entities use entity substitution with ontology XML
 output
commands:
 help print help for command
 annotate annotate ontology
 collapse minimize an ontology based on a threshold
 convert convert ontology
 diff find the differences between two ontologies
 expand expand ontology
 explain explain derivation of an inferred axiom
 export export ontology as a table
 export-prefixes export prefixes to a file
 extract extract terms from an ontology
 filter filter ontology axioms
 materialize materialize ontology
 measure compute the metrics of an ontology
 merge merge ontologies
 mirror mirror ontology imports closure
 python start a server to run ROBOT with Py4J
 query query an ontology
 reason reason ontology
 reduce reduce ontology
 relax relax ontology
 remove remove axioms from an ontology
 rename rename entities based on given mappings
 repair repair terms from an ontology
 report report terms from an ontology
 template build an ontology from a template
 unmerge unmerge ontologies
 validate-profile validate ontology against an OWL profile
 verify verify an ontology does not violate rules (as queries)
[user@cn0911 ~]$**robot -h**
usage: robot [command] [options] 
 --add-prefix  add prefix 'foo: http://bar' to the output
 --add-prefixes  add JSON-LD prefixes to the output
 --catalog  use catalog from provided file
 -h,--help print usage information
 --noprefixes do not use default prefixes
 -p,--prefix  add a prefix 'foo: http://bar'
 -P,--prefixes  use prefixes from JSON-LD file
 --strict use strict parsing when loading an ontology
 -V,--version print version information
 -v,--verbose increased logging
 -vv,--very-verbose high logging
 -vvv,--very-very-verbose maximum logging, including stack traces
 -x,--xml-entities use entity substitution with ontology XML
 output

```
 
End the interactive session:

```

[user@cn0911 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





