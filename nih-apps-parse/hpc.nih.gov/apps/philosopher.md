

document.querySelector('title').textContent = 'philosopher on Biowulf';
philosopher on Biowulf


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



Philosopher is fast, easy-to-use, scalable, and versatile data analysis software for mass spectrometry-based proteomics. Philosopher is dependency-free and can analyze both traditional database searches and open searches for post-translational modification (PTM) discovery.



### References:


* [da Veiga Leprevost, Felipe, et al. "Philosopher: a versatile toolkit for shotgun proteomics data analysis." *Nature methods* 17.9 (2020): 869-870.](https://www.nature.com/articles/s41592-020-0912-y)


Documentation
* [philosopher Main Site](https://philosopher.nesvilab.org/)
* [philosopher Documentation](https://github.com/Nesvilab/philosopher/wiki)


Important Notes
* Module Name: philosopher (see [the modules page](/apps/modules.html) for more information)
 * This application can be executed at the command line or invoked through the [fragpipe](fragpipe.html) GUI frontend.
* Raw files are processed via an extension that produces remote procedure calls (RPCs). The firewall on Biowulf compute nodes prevents this for working correctly. Files should therefore be [converted to mzML format](http://msfragger.nesvilab.org/tutorial_convert.html).



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. This example is adapted from the [tutorial here](https://github.com/Nesvilab/philosopher/wiki/Simple-Data-Analysis).   
Sample session (user input in **bold**):



```

user@biowulf ~]$ **sinteractive -c4 --mem=8g --gres=lscratch:10**
salloc.exe: Pending job allocation 5400024
salloc.exe: job 5400024 queued and waiting for resources
salloc.exe: job 5400024 has been allocated resources
salloc.exe: Granted job allocation 5400024
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0854 are ready for job
srun: error: x11: no local DISPLAY defined, skipping

[user@cn0854 ~]$ **mkdir -p /data/${USER}/philosopher\_workdir**

[user@cn0854 ~]$ **cd /data/${USER}/philosopher\_workdir**

[user@cn0854 philosopher_workdir]$ **module load fragpipe** # <-- this will load msfragger, philosopher, and add a few more things to the PATH (like mono)
[+] Loading fragpipe  14.0  on cn0854
[+] Loading msfragger  3.1.1  on cn0854
[+] Loading java 12.0.1  ...
[+] Loading philosopher  3.3.12  on cn0854
[+] Loading mono  5.18.1

[user@cn0854 philosopher_workdir]$ **cp $MSFRAGGER\_HOME/TESTDATA/\* .**

[user@cn0854 philosopher_workdir]$ **philosopher workspace --temp /lscratch/$SLURM\_JOB\_ID --init**
INFO[16:42:43] Executing Workspace  v3.3.12
INFO[16:42:43] Creating workspace
INFO[16:42:43] Done

[user@cn0854 philosopher_workdir]$ **philosopher database --id UP000005640 --contam**
INFO[16:42:54] Executing Database  v3.3.12
INFO[16:42:54] Fetching database UP000005640
INFO[16:43:09] Processing decoys
INFO[16:43:12] Creating file
INFO[16:43:24] Processing decoys
INFO[16:43:27] Creating file
INFO[16:43:28] Done

[user@cn0854 philosopher_workdir]$ **java -jar $MSFRAGGER\_JAR --config**
Creating configuration files
Writing file: /gpfs/gsfs11/users/user/philosopher_workdir/closed_fragger.params
Writing file: /gpfs/gsfs11/users/user/philosopher_workdir/open_fragger.params
Writing file: /gpfs/gsfs11/users/user/philosopher_workdir/nonspecific_fragger.params
Writing file: /gpfs/gsfs11/users/user/philosopher_workdir/Nglyco-HCD_fragger.params

[user@cn0854 philosopher_workdir]$ **ls -la**
total 1846659
drwxr-xr-x  3 user user           4096 Jan  4 16:43 .
drwxr-x--- 17 user user           4096 Jan  4 16:41 ..
-rw-r--r--  1 user user     1821252699 Jan  4 16:42 20190202_QExHFX2_RSLC8_PST_HeLa_10ng_1ulLoop_muPAC_1hr_15k_7.mzML
-rw-r--r--  1 user user       68982756 Jan  4 16:43 2021-01-04-decoys-contam-UP000005640.fas
-rw-r--r--  1 user user           8606 Jan  4 16:43 closed_fragger.params
drwxr-xr-x  2 user user           4096 Jan  4 16:43 .meta
-rw-r--r--  1 user user          12613 Jan  4 16:43 Nglyco-HCD_fragger.params
-rw-r--r--  1 user user           8601 Jan  4 16:43 nonspecific_fragger.params
-rw-r--r--  1 user user           8607 Jan  4 16:43 open_fragger.params

```


At this point you will have several configuration files in your current working directory that will determine the behavior of msfragger. Use the following sed commands to set values apporpriately. You will need to change the num\_threads argument to match the number of cpus you requested with the sinterative command above. You will also need to change the database\_name to the correct value.


```

[user@cn0854 philosopher_workdir]$ **sed -i 's/num\_threads = 0/num\_threads = 4/g' closed\_fragger.params**

[user@cn0854 philosopher_workdir]$ **sed -i 's/database\_name = test.fasta/database\_name = 2021-01-04-decoys-contam-UP000005640.fas/g' closed\_fragger.params**

[user@cn0854 philosopher_workdir]$ **sed -i 's/calibrate\_mass = 2/calibrate\_mass = 0/g' closed\_fragger.params**

```


Now you can use the edited configuration file to run msfragger and finish running philosopher. 


```

[user@cn0854 philosopher_workdir]$ **java -Xmx${SLURM\_MEM\_PER\_NODE}m -jar $MSFRAGGER\_JAR closed\_fragger.params \
 20190202\_QExHFX2\_RSLC8\_PST\_HeLa\_10ng\_1ulLoop\_muPAC\_1hr\_15k\_7.mzML**
MSFragger version MSFragger-3.1.1
Batmass-IO version 1.19.5
timsdata library version timsdata-2-7-0
(c) University of Michigan
RawFileReader reading tool. Copyright (c) 2016 by Thermo Fisher Scientific, Inc. All rights reserved.
System OS: Linux, Architecture: amd64
Java Info: 12.0.1, Java HotSpot(TM) 64-Bit Server VM, Oracle Corporation
JVM started with 8 GB memory
Checking database...
Parameter 'diagnostic_intensity_filter' was not supplied. Using default value: 0.000000
Parameter 'ion_series_definitions' was not supplied. Using default value:
Parameter 'Y_type_masses' was not supplied. Using default value:
Parameter 'diagnostic_fragments' was not supplied. Using default value:
Checking /gpfs/gsfs11/users/user/philosopher_workdir/20190202_QExHFX2_RSLC8_PST_HeLa_10ng_1ulLoop_muPAC_1hr_15k_7.mzML...

************************************MAIN SEARCH************************************
Checking database...
[snip...]
Operating on slice 2 of 2:
        Fragment index slice generated in 4.68 s
        001. 20190202_QExHFX2_RSLC8_PST_HeLa_10ng_1ulLoop_muPAC_1hr_15k_7.mzML 7.9 s | deisotoping 0.3 s
                [progress: 24044/24044 (100%) - 18481 spectra/s] 1.3s | postprocessing 10.5 s
***************************MAIN SEARCH DONE IN 1.144 MIN***************************

*******************************TOTAL TIME 1.357 MIN********************************

[user@cn0854 philosopher_workdir]$ **philosopher peptideprophet --database 2021-01-04-decoys-contam-UP000005640.fas \
 --ppm --accmass --expectscore --decoyprobs \
 --nonparam 20190202\_QExHFX2\_RSLC8\_PST\_HeLa\_10ng\_1ulLoop\_muPAC\_1hr\_15k\_7.pepXML**
INFO[16:53:41] Executing PeptideProphet  v3.3.12
 file 1: /data/user/philosopher_workdir/20190202_QExHFX2_RSLC8_PST_HeLa_10ng_1ulLoop_muPAC_1hr_15k_7.pepXML
WARNING: Unable to open file: /data/user/philosopher_workdir/20190202_QExHFX2_RSLC8_PST_HeLa_10ng_1ulLoop_muPAC_1hr_15k_7.mzML cannot correct scan numbers!
WARNING: Unable to open file: /data/user/philosopher_workdir/20190202_QExHFX2_RSLC8_PST_HeLa_10ng_1ulLoop_muPAC_1hr_15k_7.mzML cannot correct scan numbers!
 processed altogether 21514 results
INFO: Results written to file: /data/user/philosopher_workdir/interact-20190202_QExHFX2_RSLC8_PST_HeLa_10ng_1ulLoop_muPAC_1hr_15k_7.pep.xml

  - /data/user/philosopher_workdir/interact-20190202_QExHFX2_RSLC8_PST_HeLa_10ng_1ulLoop_muPAC_1hr_15k_7.pep.xml
  - Building Commentz-Walter keyword tree...
  - Searching the tree...
  - Linking duplicate entries...
  - Printing results...

using Accurate Mass Bins
using PPM mass difference
Using Decoy Label "rev_".
Decoy Probabilities will be reported.
Using non-parametric distributions
 (X! Tandem) (using Tandem's expectation score for modeling)
adding ACCMASS mixture distribution
using search_offsets in ACCMASS mixture distr: 0
init with X! Tandem stricttrypsin
MS Instrument info: Manufacturer: UNKNOWN, Model: UNKNOWN, Ionization: UNKNOWN, Analyzer: UNKNOWN, Detector: UNKNOWN

INFO: Processing standard MixtureModel ...
 PeptideProphet  (TPP v5.2.1-dev Flammagenitus, Build 201906251008-exported (Linux-x86_64)) AKeller@ISB
 read in 0 1+, 15055 2+, 5899 3+, 523 4+, 31 5+, 6 6+, and 0 7+ spectra.
Initialising statistical models ...
Found 2621 Decoys, and 18893 Non-Decoys
Iterations: .........10.........20.....
WARNING: Mixture model quality test failed for charge (1+).
WARNING: Mixture model quality test failed for charge (6+).
WARNING: Mixture model quality test failed for charge (7+).
model complete after 26 iterations
ERRO[16:54:42] unlinkat /lscratch/5400024: permission denied
INFO[16:54:42] Done

[user@cn0854 philosopher_workdir]$ **philosopher proteinprophet \
 interact-20190202\_QExHFX2\_RSLC8\_PST\_HeLa\_10ng\_1ulLoop\_muPAC\_1hr\_15k\_7.pep.xml**
INFO[16:57:02] Executing ProteinProphet  v3.3.12
ProteinProphet (C++) by Insilicos LLC and LabKey Software, after the original Perl by A. Keller (TPP v5.2.1-dev Flammagenitus, Build 201906251008-exported (Linux-x86_64))
 (no FPKM) (using degen pep info)
Reading in /data/user/philosopher_workdir/interact-20190202_QExHFX2_RSLC8_PST_HeLa_10ng_1ulLoop_muPAC_1hr_15k_7.pep.xml...
...read in 0 1+, 10738 2+, 4312 3+, 346 4+, 0 5+, 0 6+, 0 7+ spectra with min prob 0.05

Initializing 14657 peptide weights: 0%...10%...20%...30%...40%...50%...60%...70%...80%...90%...100%
Calculating protein lengths and molecular weights from database /data/user/philosopher_workdir/2021-01-04-decoys-contam-UP000005640.fas
.........:.........:.........:.........:.........:.........:.........:.........:.........:.........1000
[snip...]
Computing probabilities for 12033 proteins.  Loop 1: 0%...20%...40%...60%...80%...100%  Loop 2: 0%...20%...40%...60%...80%...100%
Computing 6023 protein groups: 0%...10%...20%...30%...40%...50%...60%...70%...80%...90%...100%
Calculating sensitivity...and error tables...
Computing MU for 12033 proteins: 0%...10%...20%...30%...40%...50%...60%...70%...80%...90%...100%
INFO: mu=6.75194e-06, db_size=104669561

Finished.
ERRO[16:58:12] unlinkat /lscratch/5400024: permission denied
INFO[16:58:12] Done

[user@cn0854 philosopher_workdir]$ **philosopher filter --razor --pepxml \
 interact-20190202\_QExHFX2\_RSLC8\_PST\_HeLa\_10ng\_1ulLoop\_muPAC\_1hr\_15k\_7.pep.xml \
 --protxml interact.prot.xml**
INFO[17:12:42] Executing Filter  v3.3.12
INFO[17:12:42] Processing peptide identification files
INFO[17:12:44] 1+ Charge profile                             decoy=0 target=0
INFO[17:12:44] 2+ Charge profile                             decoy=143 target=10595
INFO[17:12:44] 3+ Charge profile                             decoy=42 target=4270
INFO[17:12:44] 4+ Charge profile                             decoy=3 target=343
INFO[17:12:44] 5+ Charge profile                             decoy=0 target=0
INFO[17:12:44] 6+ Charge profile                             decoy=0 target=0
INFO[17:12:44] Database search results                       ions=14649 peptides=13455 psms=15396
INFO[17:12:44] Converged to 1.00 % FDR with 15147 PSMs       decoy=152 threshold=0.1032 total=15299
INFO[17:12:45] Converged to 1.00 % FDR with 13203 Peptides   decoy=133 threshold=0.1675 total=13336
INFO[17:12:45] Converged to 1.00 % FDR with 14402 Ions       decoy=145 threshold=0.1313 total=14547
INFO[17:12:48] Protein inference results                     decoy=270 target=5753
INFO[17:12:48] Converged to 1.02 % FDR with 1861 Proteins    decoy=19 threshold=0.983 total=1880
INFO[17:12:49] 2D FDR estimation: Protein mirror image       decoy=1861 target=1861
INFO[17:12:49] Second filtering results                      ions=14175 peptides=12983 psms=14916
INFO[17:12:49] Converged to 0.12 % FDR with 14898 PSMs       decoy=18 threshold=0.0513 total=14916
INFO[17:12:49] Converged to 0.13 % FDR with 12965 Peptides   decoy=18 threshold=0.053 total=12983
INFO[17:12:49] Converged to 0.12 % FDR with 14157 Ions       decoy=18 threshold=0.0519 total=14175
INFO[17:12:51] Post processing identifications
INFO[17:12:53] Mapping modifications
INFO[17:12:58] Processing protein inference
INFO[17:13:16] Assigning protein identifications to layers
INFO[17:13:17] Updating razor PSM assignment to proteins
INFO[17:13:18] Calculating spectral counts
INFO[17:13:18] Saving
ERRO[17:13:20] unlinkat /lscratch/5400024: permission denied
INFO[17:13:20] Done

[user@cn0854 philosopher_workdir]$ **philosopher filter --razor --pepxml \
 interact-20190202\_QExHFX2\_RSLC8\_PST\_HeLa\_10ng\_1ulLoop\_muPAC\_1hr\_15k\_7.pep.xml \
 --protxml interact.prot.xml**
INFO[17:12:42] Executing Filter  v3.3.12
INFO[17:12:42] Processing peptide identification files
INFO[17:12:44] 1+ Charge profile                             decoy=0 target=0
INFO[17:12:44] 2+ Charge profile                             decoy=143 target=10595
INFO[17:12:44] 3+ Charge profile                             decoy=42 target=4270
INFO[17:12:44] 4+ Charge profile                             decoy=3 target=343
INFO[17:12:44] 5+ Charge profile                             decoy=0 target=0
INFO[17:12:44] 6+ Charge profile                             decoy=0 target=0
INFO[17:12:44] Database search results                       ions=14649 peptides=13455 psms=15396
INFO[17:12:44] Converged to 1.00 % FDR with 15147 PSMs       decoy=152 threshold=0.1032 total=15299
INFO[17:12:45] Converged to 1.00 % FDR with 13203 Peptides   decoy=133 threshold=0.1675 total=13336
INFO[17:12:45] Converged to 1.00 % FDR with 14402 Ions       decoy=145 threshold=0.1313 total=14547
INFO[17:12:48] Protein inference results                     decoy=270 target=5753
INFO[17:12:48] Converged to 1.02 % FDR with 1861 Proteins    decoy=19 threshold=0.983 total=1880
INFO[17:12:49] 2D FDR estimation: Protein mirror image       decoy=1861 target=1861
INFO[17:12:49] Second filtering results                      ions=14175 peptides=12983 psms=14916
INFO[17:12:49] Converged to 0.12 % FDR with 14898 PSMs       decoy=18 threshold=0.0513 total=14916
INFO[17:12:49] Converged to 0.13 % FDR with 12965 Peptides   decoy=18 threshold=0.053 total=12983
INFO[17:12:49] Converged to 0.12 % FDR with 14157 Ions       decoy=18 threshold=0.0519 total=14175
INFO[17:12:51] Post processing identifications
INFO[17:12:53] Mapping modifications
INFO[17:12:58] Processing protein inference
INFO[17:13:16] Assigning protein identifications to layers
INFO[17:13:17] Updating razor PSM assignment to proteins
INFO[17:13:18] Calculating spectral counts
INFO[17:13:18] Saving
ERRO[17:13:20] unlinkat /lscratch/5400024: permission denied
INFO[17:13:20] Done

[user@cn0854 philosopher_workdir]$ **philosopher freequant --dir .**
INFO[17:14:35] Executing Label-free quantification  v3.3.12
INFO[17:14:38] Indexing PSM information
INFO[17:14:38] Reading spectra and tracing peaks
INFO[17:14:38] Processing 20190202_QExHFX2_RSLC8_PST_HeLa_10ng_1ulLoop_muPAC_1hr_15k_7
INFO[17:16:40] Assigning intensities to data layers
ERRO[17:16:42] unlinkat /lscratch/5400024: permission denied
INFO[17:16:42] Done

[user@cn0854 philosopher_workdir]$ **philosopher workspace --backup**
INFO[17:20:06] Executing Workspace  v3.3.12
INFO[17:20:06] Creating backup
ERRO[17:20:06] Cannot locate meta directory.
INFO[17:20:17] Done

[user@cn0854 philosopher_workdir]$ **ls ./-2021-01-04T17\:20\:06-05\:00.zip**
./-2021-01-04T17:20:06-05:00.zip

[user@cn0854 philosopher_workdir]$ **exit**
exit
salloc.exe: Relinquishing job allocation 5400024
salloc.exe: Job allocation 5400024 has been revoked.

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. philosopher.sh). For example:



```

#!/bin/bash
set -e
module load fragpipe
cd /data/${USER}/workdir
java -Xmx${SLURM_MEM_PER_NODE}m -jar $MSFRAGGER_JAR closed_fragger.params \
    20190202_QExHFX2_RSLC8_PST_HeLa_10ng_1ulLoop_muPAC_1hr_15k_7.mzML
philosopher peptideprophet --database 2021-01-04-decoys-contam-UP000005640.fas \
    --ppm --accmass --expectscore --decoyprobs \
    --nonparam 20190202_QExHFX2_RSLC8_PST_HeLa_10ng_1ulLoop_muPAC_1hr_15k_7.pepXML
philosopher proteinprophet \
    interact-20190202_QExHFX2_RSLC8_PST_HeLa_10ng_1ulLoop_muPAC_1hr_15k_7.pep.xml
philosopher filter --razor --pepxml \
    interact-20190202_QExHFX2_RSLC8_PST_HeLa_10ng_1ulLoop_muPAC_1hr_15k_7.pep.xml \
    --protxml interact.prot.xml
philosopher freequant --dir .

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] philosopher.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. philosopher.swarm). For example:



```

philosopher peptideprophet --database 2021-01-04-decoys-contam-UP000005640.fas \
    --ppm --accmass --expectscore --decoyprobs \
    --nonparam file1.pepXML
philosopher peptideprophet --database 2021-01-04-decoys-contam-UP000005640.fas \
    --ppm --accmass --expectscore --decoyprobs \
    --nonparam file2.pepXML
philosopher peptideprophet --database 2021-01-04-decoys-contam-UP000005640.fas \
    --ppm --accmass --expectscore --decoyprobs \
    --nonparam file3.pepXML

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f philosopher.swarm [-g #] [-t #] --module fragpipe
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module fragpipe Loads the fragpipe module for each subjob in the swarm 
 | |
 | |
 | |








