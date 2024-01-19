

document.querySelector('title').textContent = "REMORA";
REMORA on Biowulf


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



Methylation/modified base calling separated from basecalling. Remora primarily provides an API to call modified bases for basecaller programs such as Bonito. Remora also provides the tools to prepare datasets, train modified base models and run simple inference.



### References:


* Please see:
 [**Original documentation**](https://github.com/nanoporetech/remora)
*Oxford Nanopore Technologies*


Documentation
* [REMORA Main Site](https://github.com/nanoporetech/remora)


Important Notes
* Module Name: remora (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded
* Environment variables set 
	+ REMORA\_HOME* Example files in /usr/local/apps/remora/2.1.1/tests/data



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
[user@cn4338 ~]$ **module load remora**
[+] Loading remora  2.1.1  on cn4338 
[+] Loading singularity  3.10.5  on cn4338

```


Example
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Running help command:



```

[user@cn4338] **cp -a /usr/local/apps/remora/2.1.1/tests/data .**
[user@cn4338 data]$ **remora --help**
remora --help
usage: remora [-h] [-v] {dataset,model,infer,validate,analyze} ...
  
********** Remora *********
  
Modified base model training and application.
  
optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Show Remora version and exit.
  
sub-commands:
  dataset             Create or perform operations on a Remora dataset
  model               Train or perform operations on Remora models
  infer               Perform Remora model inference
  validate            Validate modified base predictions
  analyze             Analyze nanopore data including raw signal

```

Generate Training Data:



```

[apptest1@cn4338 data]$ **remora \
 dataset prepare \
 can\_reads.pod5 \
 can\_mappings.bam \
 --output-remora-training-file can\_chunks.npz \
 --log-filename prep\_can.log \
 --refine-kmer-level-table levels.txt \
 --refine-rough-rescale \
 --motif CG 0 \
 --mod-base-control** 
Indexing BAM by read id: 14 Reads [00:00, 10146.93 Reads/s]
[14:30:14] Extracting read IDs from POD5
[14:30:14] Found 14 BAM records, 14 POD5 reads, and 14 in common
[14:30:14] Making reference-anchored training data
[14:30:14] Allocating memory for output tensors
[14:30:14] Processing reads
Extracting chunks: 100%|█████████████████████████| 14/14 [00:00<00:00, 163.20 Reads/s]
[14:30:14] Extracted 205 chunks from 14 reads.
[14:30:14] Label distribution: Counter({0: 205})

```



|  |  |
| --- | --- |
|  For more examples, please visit the [Remora Github Page](https://github.com/nanoporetech/remora) | |








