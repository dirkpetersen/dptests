

document.querySelector('title').textContent = "TEMPLATE";
NESTED on Biowulf


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



nested (now also called TE-greedy) is software to analyze nested LTR transposable elements
in DNA sequences, such as reference genomes. It is made of two components: nested-generator
for generating simulated sequences of nested retrotransposons, and nested-nester (now called
TE-greedy-nester) that looks for nested, as well as non-nested and solo-LTR repeat sequences
in the input. Unlike other similar software, TE-greedy-nester is structure-based by using
de-novo retrotransposon identification software LTR Finder, relying on sequence information
only secondarily.



### References:


* Matej Lexa, Pavel Jedlicka, Ivan Vanat, Michal Cervenansky, Eduard Kejnovsky
 [**TE-greedy-nester: structure-based detection of LTR retrotransposons and their nesting**](https://doi.org/10.1093/bioinformatics/btaa632)
*Bioinformatics, Volume 36, Issue 20, October 2020, Pages 4991–4999*


Documentation
* [Chap Main Site](https://gitlab.fi.muni.cz/lexa/nested)


Important Notes
* Module Name: nested (see [the modules page](/apps/modules.html) for more information)

* Example files in /usr/local/apps/nested/2.0.0/test\_data



Getting Started
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

[user@cn3144 ~]$ **module load nested**
[+] Loading nested  2.0.0  on cn3144
[+] Loading singularity  3.10.5  on cn3144

[user@cn3144 ~]$**nested-nester --help**
Usage: nested-nester [OPTIONS] INPUT_FASTA

Options:
  -s, --sketch                    Sketch output.
  -f, --format TEXT               Format for GFF.
  -o, --output_fasta_offset INTEGER
                                  Number of bases around the element included
                                  in output fasta files.
  -d, --output_folder PATH        Output data folder.
  -t, --initial_threshold INTEGER
                                  Initial threshold value.
  -m, --threshold_multiplier FLOAT
                                  Threshold multiplier.
  -n, --threads INTEGER           Number of threads
  -dt, --discovery_tool [LTR_finder|LTRharvest|finder|harvest]
                                  Determines which tool is used for
                                  retrotransoson discovery. Default:
                                  LTR_finder
  -solo, --solo_ltrs              Run solo LTR module
  --help                          Show this message and exit.

[user@cn3144 ~]$ **nested-generator --help**
Usage: nested-generator [OPTIONS] INPUT_DB OUTPUT_DB

Options:
  -l, --baselength INTEGER        Baselength for generated elements.
  -i, --number_of_iterations INTEGER
                                  Number of inserted elements.
  -n, --number_of_elements INTEGER
                                  Number of generated sequences.
  -f, --filter                    Filter database and create new one with
                                  given output db path.
  -s, --filter_string TEXT        Filter entries by given string [ONLY
                                  RELEVANT WITH -filter OPTION].
  -o, --filter_offset INTEGER     LTR offset allowed [ONLY RELEVANT WITH
                                  -filter OPTION].
  -p, --percentage INTEGER        Percentage of elements in generated
                                  sequence.
  -a, --average_element INTEGER   Average element length in database.
  -e, --expected_length INTEGER   Expected output sequence length [ONLY WORKS
                                  WITH -percentage and -average_element].
  -d, --output_directory TEXT     Output directory.
  --help                          Show this message and exit.

```


Example
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Example running nested



```

[user@cn3144 ~]$**cp -a /usr/local/apps/nested/2.0.0/test\_data/ .** 
[user@cn3144 ~]$ **nested-nester test\_data/151kb\_adh1\_bothPrimes.fasta**
/usr/local/lib/python3.9/dist-packages/nested-1.0.0-py3.9.egg/nested/config/config.py:14: YAMLLoadWarning: calling yaml.load() without Loader=... is deprecated, as the default Loader is unsafe. Please read https://msg.pyyaml.org/load for full details.
Processing adh1_vicinity_150bpPlMin
Processing adh1_vicinity_150bpPlMin: DONE [0:00:10.348563]
Total time: 0:00:10.350736
Number of errors: 0

```

For more information please see the [GitLab Page](https://gitlab.fi.muni.cz/lexa/nested/-/blob/master/README.md) 








