

document.querySelector('title').textContent = 'Cooltools: enabling high-resolution Hi-C analysis in Python';
**Cooltools: enabling high-resolution Hi-C analysis in Python**


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



Cooltools is a suite of computational tools that enables flexible, scalable, and
reproducible analysis of high-resolution contact frequency data. Cooltools leverages the widely-adopted
cooler format which handles storage and access for high-resolution datasets. Cooltools provides a paired
command line interface and Python application programming interface, which respectively
facilitate workflows on high-performance computing clusters and in interactive analysis environments.



### References:


* Open2C, Nezar Abdennur, Sameer Abraham, Geoffrey Fudenberg, Ilya M. Flyamer, Aleksandra A.
Galitsyna, Anton Goloborodko, Maxim Imakaev, Betul A. Oksuz, Sergey V. Venev,   

*Cooltools: enabling high-resolution Hi-C analysis in Python*   

[Nature Communications](https://www.nature.com/articles/s41467-021-24041-8) **12**, Article number: 3836 (2021).


Documentation
* [cooltools Github page](https://github.com/open2c/cooltools)
* [cooltools documentation](https://cooltools.readthedocs.io/en/latest/#)


Important Notes
* Module Name: cooltools (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **CT\_HOME**  installation directory
	+ **CT\_BIN**       executable directory
	+ **CT\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive** 
[user@cig 3335 ~]$ **module load cooltools**
[+] Loading singularity  3.10.5  on cn3335
[+] Loading cooltools  0.5.4
[user@cn3335 ~]$ **cooltools**
Usage: cooltools [OPTIONS] COMMAND [ARGS]...

  Type -h or --help after any subcommand for more information.

Options:
  -v, --verbose  Verbose logging
  -d, --debug    Post mortem debugging
  -V, --version  Show the version and exit.
  -h, --help     Show this message and exit.

Commands:
  coverage        Calculate the sums of cis and genome-wide contacts (aka...
  dots            Call dots on a Hi-C heatmap that are not larger than...
  eigs-cis        Perform eigen value decomposition on a cooler matrix to...
  eigs-trans      Perform eigen value decomposition on a cooler matrix to...
  expected-cis    Calculate expected Hi-C signal for cis regions of...
  expected-trans  Calculate expected Hi-C signal for trans regions of...
  genome          Utilities for binned genome assemblies.
  insulation      Calculate the diamond insulation scores and call...
  pileup          Perform retrieval of the snippets from .cool file.
  random-sample   Pick a random sample of contacts from a Hi-C map.
  saddle          Calculate saddle statistics and generate saddle plots...
  virtual4c       Generate virtual 4C profile from a contact map by...
[user@cn3335 ~]$ **ct\_python**
Python 3.8.18 (default, Sep 11 2023, 13:40:15)
[GCC 11.2.0] :: Anaconda, Inc. on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import cooltools
>>>

```





