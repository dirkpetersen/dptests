

document.querySelector('title').textContent = 'Autodock Crankprep';
Autodock Crankprep


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



 AutoDock CrankPep or ADCP is an AutoDock docking engine specialized for docking peptides. 
 It combines technology form the protein folding filed with an efficient representation of
 a rigid receptor as affinity grids to fold the peptide in the context of the energy landscape
 created by the receptor.



### References:


* Zhang, Y.; Sanner, M. F.
 [**AutoDock CrankPep: Combining Folding and Docking to Predict Protein-Peptide Complexes.**](https://doi.org/10.1093/bioinformatics/btz459)
*Bioinformatics, Volume 35, Issue 24, December 2019, Pages 5121–5127*


Documentation
* [AutodockCrankPrep Main Site](https://ccsb.scripps.edu/adcp/)


Important Notes
* Module Name: AutodockCrankprep (see [the modules page](/apps/modules.html) for more information)
* Example files in /usr/local/apps/AutodockCrankprep/data



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

[user@cn3144 ~]$ **module load AutodockCrankprep**
[+] Loading AutodockCrankprep  1.0  on cn2882
[user@cn3144 ~]$ **prepare\_receptor**
prepare_receptor4: receptor filename must be specified.
Usage: prepare_receptor4.py -r filename

  Description of command...
    -r  receptor_filename 
    supported file types include pdb,mol2,pdbq,pdbqs,pdbqt, possibly pqr,cif
  Optional parameters:
    [-v] verbose output (default is minimal output)
    [-o pdbqt_filename]  (default is 'molecule_name.pdbqt')
    [-A] type(s) of repairs to make: 
      'bonds_hydrogens': build bonds and add hydrogens 
      'bonds': build a single bond from each atom with no bonds to its closest neighbor
      'hydrogens': add hydrogens
      'checkhydrogens': add hydrogens only if there are none already
      'None': do not make any repairs 
      (default is 'None')
    [-C] preserve all input charges ie do not add new charges 
      (default is addition of gasteiger charges)
    [-p] preserve input charges on specific atom types, eg -p Zn -p Fe
    [-U] cleanup type:
      'nphs': merge charges and remove non-polar hydrogens
      'lps': merge charges and remove lone pairs
      'waters': remove water residues
      'nonstdres': remove chains composed entirely of residues of
              types other than the standard 20 amino acids
      'deleteAltB': remove XX@B atoms and rename XX@A atoms->XX
              (default is 'nphs_lps_waters_nonstdres') 
    [-e] delete every nonstd residue from any chain
      'True': any residue whose name is not in this list:
              ['CYS','ILE','SER','VAL','GLN','LYS','ASN', 
              'PRO','THR','PHE','ALA','HIS','GLY','ASP', 
              'LEU', 'ARG', 'TRP', 'GLU', 'TYR','MET', 
              'HID', 'HSP', 'HIE', 'HIP', 'CYX', 'CSS']
      will be deleted from any chain. 
      NB: there are no  nucleic acid residue names at all 
      in the list and no metals. 
      (default is False which means not to do this)
    [-M]  interactive 
            (default is 'automatic': outputfile is written with no further user input)
    [-d dictionary_filename] file to contain receptor summary information
    [-w]   assign each receptor atom a unique name: newname is original name plus its index(1-based)

```


Example
Create a batch input file (e.g. TEMPLATE.sh). For example:



```

#!/bin/bash
set -e
module load TEMPLATE
TEMPLATE < TEMPLATE.in > TEMPLATE.out

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] TEMPLATE.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. TEMPLATE.swarm). For example:



```

TEMPLATE < TEMPLATE.in > TEMPLATE.out
TEMPLATE < TEMPLATE.in > TEMPLATE.out
TEMPLATE < TEMPLATE.in > TEMPLATE.out
TEMPLATE < TEMPLATE.in > TEMPLATE.out

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f TEMPLATE.swarm [-g #] [-t #] --module TEMPLATE
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module TEMPLATE Loads the TEMPLATE module for each subjob in the swarm 
 | |
 | |
 | |








