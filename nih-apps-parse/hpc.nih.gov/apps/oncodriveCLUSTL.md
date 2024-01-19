

document.querySelector('title').textContent = 'OncodriveCLUSTL: a sequence-based clustering method to identify cancer drivers';
**OncodriveCLUSTL: a sequence-based clustering method to identify cancer drivers**


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



OncodriveCLUSTL is a sequence-based clustering algorithm to detect significant clustering signals 
across genomic regions. 
It is based based on a local background model 
derived from the simulation of mutations 
accounting for the composition of trior penta-nucleotide context substitutions observed in the cohort under study.



### References:


* Claudia Arnedo-Pac, Loris Mularoni, Ferran Muinos, Abel Gonzalez-Perez and Nuria Lopez-Bigas,   

*OncodriveCLUSTL: a sequence-based clustering method to identify cancer drivers*   

[Bioinformatics](https://academic.oup.com/bioinformatics/article/35/22/4788/5522012), 2019 **35**(22), 4788–4790. doi: 10.1093/bioinformatics/btz501


Documentation
* [OncodriveCLUSTL Bitbucket repository](https://bitbucket.org/bbglab/oncodriveclustl/src/master/)
* [OncodriveCLUSTL Documentation](https://pypi.org/project/oncodriveclustl/)


Important Notes
* Module Name: oncodriveCLUSTL (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **ODCLUSTL\_HOME**  installation directory
	+ **ODCLUSTL\_BIN**       executable directory
	+ **ODCLUSTL\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=75g -c20 --gres=lscratch:20**
[user@cn3335 ~]$ **module load oncodriveCLUSTL**
[+] Loading singularity  3.10.5  on cn4338
[+] Loading oncodriveCLUSTL  1.1.1
[user@cn3335 ~]$ **oncodriveclustl -h**
Usage: oncodriveclustl [OPTIONS]

  OncodriveCLUSTL is a sequence based clustering method to identify cancer
  drivers across the genome

  Args:     input_file (str): path to mutations file     regions_file (str):
  path to input genomic coordinates file     output_directory(str): path to
  output directory. Output files will be generated in it.
  input_signature (str): path to file containing input context based
  mutational probabilities.         By default (when no input signatures),
  OncodriveCLUSTL will calculate them from the mutations input file.
  elements_file (str): path to file containing one element per row
  (optional) to analyzed the listed elements.         By default,
  OncodriveCLUSTL analyzes all genomic elements contained in `regions_file`.
  elements (str): genomic element symbol (optional). The analysis will be
  performed only on the specified GEs.     genome (str): genome to use:
  'hg38', 'hg19', 'mm10', 'c3h', 'car', 'cast' and 'f344'
  element_mutations (int): minimum number of mutations per genomic element
  to undertake analysis     cluster_mutations (int): minimum number of
  mutations to define a cluster     smooth_window (int): Tukey kernel
  smoothing window length     cluster_window (int): clustering window length
  kmer (int): context nucleotides to calculate the mutational probabilities
  (trinucleotides or pentanucleotides)     n_simulations (int): number of
  simulations     simulation_mode (str): simulation mode
  simulation_window (int): window length to simulate mutations
  signature_calculation (str): signature calculation, mutation frequencies
  (default) or mutation counts         normalized by k-mer region counts
  signature_group (str): header of the column to group signatures. One
  signature will be computed for each group     cores (int): number of CPUs
  to use     seed (int): seed     log_level (str): verbosity of the logger
  concatenate (bool): flag to calculate clustering on collapsed genomic
  regions (e.g., coding regions in a gene)     clustplot (bool): flag to
  generate a needle plot with clusters for an element     qqplot (bool):
  flat to generate a quantile-quantile (QQ) plot for a dataset     gzip
  (bool): flag to generate GZIP compressed output files

  Returns:     None

Options:
  -i, --input-file PATH           File containing somatic mutations
                                  [required]
  -r, --regions-file PATH         File with the genomic regions to analyze
                                  [required]
  -o, --output-directory TEXT     Output directory to be created  [required]
  -sig, --input-signature PATH    File containing input context based
                                  mutational probabilities (signature)
  -ef, --elements-file PATH       File with the symbols of the elements to
                                  analyze
  -e, --elements TEXT             Symbol of the element(s) to analyze
  -g, --genome [hg38|hg19|mm10|c3h|car|cast|f344]
                                  Genome to use
  -emut, --element-mutations INTEGER
                                  Cutoff of element mutations. Default is 2
  -cmut, --cluster-mutations INTEGER
                                  Cutoff of cluster mutations. Default is 2
  -sw, --smooth-window INTEGER RANGE
                                  Smoothing window. Default is 11
  -cw, --cluster-window INTEGER RANGE
                                  Cluster window. Default is 11
  -kmer, --kmer [3|5]             K-mer nucleotide context
  -n, --n-simulations INTEGER     number of simulations. Default is 1000
  -sim, --simulation-mode [mutation_centered|region_restricted]
                                  Simulation mode
  -simw, --simulation-window INTEGER RANGE
                                  Simulation window. Default is 31
  -sigcalc, --signature-calculation [frequencies|region_normalized]
                                  Signature calculation: mutation frequencies
                                  (default) or k-mer mutation counts
                                  normalized by k-mer region counts
  -siggroup, --signature-group [SIGNATURE|SAMPLE|CANCER_TYPE]
                                  Header of the column to group signatures
                                  calculation
  -c, --cores INTEGER RANGE       Number of cores to use in the computation.
                                  By default it will use all the available
                                  cores.
  --seed INTEGER                  Seed to use in the simulations
  --log-level [debug|info|warning|error|critical]
                                  Verbosity of the logger
  --concatenate                   Calculate clustering on concatenated genomic
                                  regions (e.g., exons in coding sequences)
  --clustplot                     Generate a needle plot with clusters for an
                                  element
  --qqplot                        Generate a quantile-quantile (QQ) plot for a
                                  dataset
  --gzip                          Gzip compress files
  -h, --help                      Show this message and exit.

```

Copy sample data to the current folder:

```

[user@cn3335 ~]$ **cp -r $ODCLUSTL\_DATA/\* .**

```

Now let's run oncodriveCLUSTL on the sample data. According to the the oncodriveCLUSTL documentation,
"The first time that you run OncodriveCLUSTL with a given reference genome, it will download it from our servers. By default the downloaded datasets go to ~/.bgdata. If you want to move these datasets to another folder you have to define the system environment variable BGDATA\_LOCAL with an export command."

```

[user@cn3335 ~]$ **oncodriveclustl -i PAAD.tsv.gz -r cds.hg19.regions.gz -o test\_output**
2023-02-02 08:32:50,073 [110140] INFO     root: OncodriveCLUSTL
2023-02-02 08:32:50,073 [110140] INFO     root:
input_file: PAAD.tsv.gz
regions_file: cds.hg19.regions.gz
input_signature: None
output_directory: test_output
genome: hg19
element_mutations: 2
cluster_mutations: 2
concatenate: False
smooth_window: 11
cluster_window: 11
k-mer: 3
simulation_mode: mutation_centered
simulation_window: 31
n_simulations: 1000
signature_calculation: frequencies
signature_group: None
cores: 128
gzip: False
seed: None
2023-02-02 08:32:50,075 [110140] INFO     root: Initializing OncodriveCLUSTL...
2023-02-02 08:32:50,077 [110140] WARNING  root:
Running with default simulating, smoothing and clustering OncodriveCLUSTL parameters. Default parameters may not be optimal for your data.
Please, read Supplementary Methods to perform model selection for your data.
2023-02-02 08:32:50,079 [110140] WARNING  root:
Signatures will be calculated as mutation frequencies: # mutated ref>alt k-mer counts / # total substitutions
Please, read Supplementary Methods to perform a more accurate signatures calculation
2023-02-02 08:32:50,080 [110140] INFO     root: Parsing genomic regions and mutations...
2023-02-02 08:33:01,448 [110140] INFO     root: Regions parsed
2023-02-02 08:33:01,639 [110140] INFO     root: Mutations parsed
2023-02-02 08:33:01,714 [110140] INFO     root: Validated elements in genomic regions: 20169
2023-02-02 08:33:01,715 [110140] INFO     root: Validated elements with mutations: 5183
2023-02-02 08:33:01,716 [110140] INFO     root: Total substitution mutations: 7913
2023-02-02 08:33:01,717 [110140] INFO     root: Computing signature...
2023-02-02 08:33:05,327 [110140] INFO     root: Signature computed
2023-02-02 08:33:05,349 [110140] INFO     root: Calculating results 1456 elements...
2023-02-02 08:33:05,352 [110140] INFO     root: Iteration 1 of 15
                   simulations: 100%|█████████████████████████████████| 3/3 [00:14<00:00,  4.84s/it]
               post processing:  79%|██████████████████████▉      | 101/128 [01:36<00:25,  1.04it/s]
2023-02-02 08:34:59,477 [110140] INFO     root: Iteration 2 of 15
                   simulations: 100%|█████████████████████████████████| 7/7 [00:13<00:00,  1.95s/it]
               post processing:  79%|██████████████████████▉      | 101/128 [01:30<00:24,  1.11it/s]
2023-02-02 08:36:45,989 [110140] INFO     root: Iteration 3 of 15
                   simulations: 100%|█████████████████████████████████| 5/5 [00:22<00:00,  4.46s/it]
               post processing:  79%|██████████████████████▉      | 101/128 [01:17<00:20,  1.30it/s]
2023-02-02 08:38:29,579 [110140] INFO     root: Iteration 4 of 15
                   simulations: 100%|█████████████████████████████████| 5/5 [00:18<00:00,  3.70s/it]
               post processing:  79%|██████████████████████▉      | 101/128 [01:34<00:25,  1.07it/s]
2023-02-02 08:40:25,999 [110140] INFO     root: Iteration 5 of 15
                   simulations: 100%|█████████████████████████████████| 7/7 [00:30<00:00,  4.41s/it]
               post processing:  79%|██████████████████████▉      | 101/128 [01:38<00:26,  1.02it/s]
2023-02-02 08:42:39,555 [110140] INFO     root: Iteration 6 of 15
                   simulations: 100%|███████████████████████████████| 12/12 [00:30<00:00,  2.54s/it]
               post processing:  79%|██████████████████████▉      | 101/128 [01:16<00:20,  1.32it/s]
2023-02-02 08:44:30,511 [110140] INFO     root: Iteration 7 of 15
                   simulations: 100%|█████████████████████████████████| 7/7 [00:15<00:00,  2.23s/it]
               post processing:  79%|██████████████████████▉      | 101/128 [01:25<00:22,  1.19it/s]
2023-02-02 08:46:15,403 [110140] INFO     root: Iteration 8 of 15
                   simulations: 100%|███████████████████████████████| 11/11 [00:17<00:00,  1.59s/it]
               post processing:  79%|██████████████████████▉      | 101/128 [01:29<00:23,  1.13it/s]
2023-02-02 08:48:07,030 [110140] INFO     root: Iteration 9 of 15
                   simulations: 100%|███████████████████████████████| 13/13 [00:51<00:00,  3.95s/it]
               post processing:  79%|██████████████████████▉      | 101/128 [01:34<00:25,  1.07it/s]
2023-02-02 08:50:36,900 [110140] INFO     root: Iteration 10 of 15
                   simulations: 100%|█████████████████████████████████| 9/9 [00:13<00:00,  1.53s/it]
               post processing:  79%|██████████████████████▉      | 101/128 [01:41<00:27,  1.00s/it]
2023-02-02 08:52:35,327 [110140] INFO     root: Iteration 11 of 15
                   simulations: 100%|█████████████████████████████████| 5/5 [00:14<00:00,  2.96s/it]
               post processing:  79%|██████████████████████▉      | 101/128 [01:28<00:23,  1.14it/s]
2023-02-02 08:54:22,365 [110140] INFO     root: Iteration 12 of 15
                   simulations: 100%|███████████████████████████████| 10/10 [00:20<00:00,  2.09s/it]
               post processing:  79%|██████████████████████▉      | 101/128 [01:37<00:25,  1.04it/s]
2023-02-02 08:56:24,457 [110140] INFO     root: Iteration 13 of 15
                   simulations: 100%|█████████████████████████████████| 7/7 [00:38<00:00,  5.56s/it]
               post processing:  79%|██████████████████████▉      | 101/128 [01:28<00:23,  1.14it/s]
2023-02-02 08:58:35,490 [110140] INFO     root: Iteration 14 of 15
                   simulations: 100%|███████████████████████████████| 10/10 [00:18<00:00,  1.85s/it]
               post processing:  79%|██████████████████████▉      | 101/128 [01:33<00:24,  1.09it/s]
2023-02-02 09:00:31,237 [110140] INFO     root: Iteration 15 of 15
                   simulations: 100%|█████████████████████████████████| 1/1 [00:05<00:00,  5.99s/it]
               post processing:  45%|█████████████▎                | 57/128 [00:56<01:10,  1.01it/s]
2023-02-02 09:01:40,325 [110140] INFO     root: Elements results calculated
2023-02-02 09:01:40,381 [110140] INFO     root: Clusters results calculated
2023-02-02 09:01:40,383 [110140] INFO     root: Finished

```





