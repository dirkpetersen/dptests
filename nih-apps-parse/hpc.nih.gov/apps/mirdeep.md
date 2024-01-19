

document.querySelector('title').textContent = 'mirdeep on Biowulf';
mirdeep on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



miRDeep2 uses the distribution of next generation sequencing reads in the genome
along with RNA structure prediction to discover and quantitate the expression
of known and novel miRNAs. miRDeep2 represents a complete overhaul of the original
miRDeep tool.




miRDeep2 is a collection of perl scripts tied together by 3 main scripts:
* `miRDeep2.pl` - Wrapper function for the miRDeep2.pl 
 program package
* `mapper.pl` - Processes reads and/or maps them to the 
 reference genome
* `quantifier.pl` - The module maps the deep sequencing 
 reads to predefined miRNA precursors and determines by that the 
 expression of the corresponding miRNAs.


Of these, `mapper.pl` and `quantifier.pl` may run 
multithreaded bowtie suprocesses. `-o` determines the 
thread count for `mapper.pl`, and `-T` for 
`quantifier.pl`. Because of this mixed nature of processes,
it's best to run individual steps separately rather than combining them
into a single batch script.



### References:


* Marc R. Friedländer, S. D. Mackowiak, N. Li, W. Chen, and N. Rajewsky.
 *miRDeep2 accurately identifies known and hundreds of novel microRNA 
 genes in seven animal clades.* Nucleic Acids Res. 2012, 40:37-52.
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/21911355)
  | 
 [PMC](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3245920/)
  | 
 [Journal](http://nar.oxfordjournals.org/content/40/1/37.long)
* Sebastian D. Mackowiak. *UNIT 12.10 Identification of Novel and 
 Known miRNAs in Deep-Sequencing Data with miRDeep2.*
 Curr Protoc Bioinformatics. 2011, Unit 12.10.
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/22161567)
  | 
 PMC
  | 
 [Journal](http://onlinelibrary.wiley.com/doi/10.1002/0471250953.bi1210s36/abstract;jsessionid=334E68086BBEDFDD226424C4DD4FA87D.f01t03)


Documentation
* [Home page](https://www.mdc-berlin.de/8551903/en/)
* [Documentation](https://www.mdc-berlin.de/36105849/en/research/research_teams/systems_biology_of_gene_regulatory_elements/projects/miRDeep/documentation)


Important Notes
* Module Name: mirdeep (see [the modules page](/apps/modules.html) for more information)
* Example files in `$MIRDEEP_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$

[user@cn3144 ~]$ **module load mirdeep**
[user@cn3144 ~]$ **cp -r $MIRDEEP\_TEST\_DATA tutorial**
[user@cn3144 ~]$ **cd tutorial**
[user@cn3144 ~]$ # the following step is not necessary if a prebuilt
[user@cn3144 ~]$ # genome index is available
[user@cn3144 ~]$ **bowtie-build cel\_cluster.fa cel\_cluster**

[user@cn3144 ~]$ **mapper.pl reads.fa -c -j -k TCGTATGCCGTCTTCTGCTTGT \
 -l 18 -m -p cel\_cluster -o 2 -v -n \
 -s reads\_collapsed.fa \
 -t reads\_collapsed\_vs\_genome.arf**
.....
[user@cn3144 ~]$ **miRDeep2.pl reads\_collapsed.fa cel\_cluster.fa \
 reads\_collapsed\_vs\_genome.arf \
 mature\_ref\_this\_species.fa \
 mature\_ref\_other\_species.fa \
 precursors\_ref\_this\_species.fa \
 -t C.elegans 2> report.log**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).

The single threaded and multi threaded steps of the miRDeep2 pipeline could
be tied together with [snakemake](/apps/snakemake.html) or a 
similar workflow tool capable of sumbitting batch jobs.
For the example here, we will simply write a script that uses job
dependencies to tie together the three steps of the whole pipeline:




```

#! /bin/bash
cd /data/$USER/test_data/mirdeep
cp -r /usr/local/apps/mirdeep/2.0.0.7/tutorial .
cd tutorial
module load mirdeep
bowtie-build cel_cluster.fa cel_cluster

# create files for each job in the pipeline
cat > step1.sh <<'EOF'
#! /bin/bash
#SBATCH --job-name=mirdeep_s1
mapper.pl reads.fa -c -j -k TCGTATGCCGTCTTCTGCTTGT  \
  -l 18 -m -p cel_cluster -v \
  -o ${SLURM_CPUS_PER_TASK} \
  -s reads_collapsed.fa \
  -t reads_collapsed_vs_genome.arf
EOF

cat > step2.sh <<'EOF'
#! /bin/bash
#SBATCH --job-name=mirdeep_s2
quantifier.pl -p precursors_ref_this_species.fa \
  -m mature_ref_this_species.fa \
  -T ${SLURM_CPUS_PER_TASK} \
  -r reads_collapsed.fa -t cel -y 16_19
EOF

cat > step3.sh <<'EOF'
#! /bin/bash
#SBATCH --job-name=mirdeep_s3
miRDeep2.pl reads_collapsed.fa cel_cluster.fa \
  reads_collapsed_vs_genome.arf \
  mature_ref_this_species.fa \
  mature_ref_other_species.fa \
  precursors_ref_this_species.fa \
  -t C.elegans 2> report.log
EOF

# set up the pipeline run
jid1=$(sbatch -c4 step1.sh)
jid2=$(sbatch -c4 --dependency=afterany:${jid1} step2.sh)
jid3=$(sbatch --dependency=afterany:${jid2} step3.sh)

```


The script above will submit all three steps as separate jobs. Each job will
only execute if the previous job finished successfully. Of course, the same
can be achieved by manually creating batch scripts for each job and sumbitting
them individually as batch jobs.









