

document.querySelector('title').textContent = 'eggNOG-mapper: fast genome-wide functional annotation through orthology assignment. ';
**eggNOG-mapper: fast genome-wide functional annotation through orthology assignment.** 


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


  


eggNOG-mapper is a tool for functional annotation of large sets of sequences based on fast orthology assignments using precomputed clusters and phylogenies from the eggNOG database. 
Orthology assignment is ideally suited for functional inference. 
However, predicting orthology is computationally intensive at large scale, and most 
other pipelines are relatively inaccessible (e.g., new assignments only available through database updates), so less precise homology-based functional transfer was previously the default for (meta-)genome annotation. 



  

### References:


* Carlos P. Cantalapiedra, Ana Hernandez-Plaza, Ivica Letunic, Peer Bork and Jaime Huerta-Cepas.   

 *eggNOG-mapper v2: Functional Annotation, Orthology Assignments, and Domain Prediction at the Metagenomic Scale.*   

 [Molecular Biology and Evolution](https://academic.oup.com/mbe/article/38/12/5825/6379734) **38i**(12):5825–5829.
* J.Huerta-Cepas, D.Szklarczyk, D.Heller, A.Hernández-Plaza, S.K.Forslund, H.Cook, 
 D.R.Mende, I.Letunic, T.Rattei, L.J.Jensen, C. von Mering, P.Bork   

 *eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource based on 5090 organisms and 2502 viruses.*   

 [Nucleic Acids Research,](https://academic.oup.com/nar/article/47/D1/D309/5173662) Volume 47, Issue D1, 08 January 2019, Pages D309–D314.


Documentation
* [eggNOG-mapper GitHub page](https://github.com/eggnogdb/eggnog-mapper)
* [eggNOG-mapper Documentation](https://github.com/eggnogdb/eggnog-mapper/wiki)


Important Notes
* Module Name: eggNOG-mapper (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **EGGNOG\_HOME**  installation directory
	+ **EGGNOG\_BIN**       executable directory
	+ **EGGNOG\_SRC**       source code directory
	+ **EGGNOG\_TEST**  sample data directory
	 + **EGGNOG\_DATA\_DIR**  database directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=16g --gres=lscratch:10**
[user@cn3200 ~]$ **module load eggnog-mapper** 
[+] Loading python 3.7  ...
[+] Loading eggnog-mapper 2.1.2  ...
[user@cn3200 ~]$**ls $EGGNOG\_BIN**
create_dbs.py  download_eggnog_data.py  hmm_mapper.py  hmm_worker.py  shell
diamond        emapper.py               hmm_server.py  python

```

Download test data files:

```

[user@cn3200 ~]$ **mkdir /data/$USER/eggnog-mapper && cd /data/$USER/eggnog-mapper**
[user@cn3200 ~]$ **cp $EGGNOG\_TEST/\* .**

```

Run emapper.py on the test data using bacteria.dmnd diamond database:

```

[user@cn3200 ~]$**emapper.py --dmnd\_db $EGGNOG\_DATA\_DIR/bacteria.dmnd -i test\_queries.fa -o test**
 /opt/conda/envs/eggNOGmapper/lib/python3.7/site-packages/eggnog_mapper-2.1.6-py3.7.egg/eggnogmapper/bin/diamond blastp -d /usr/local/apps/eggnog-mapper/2.1.6/data/bacteria.dmnd -q /gs7/users/user/eggnog-mapper/test_queries.fa --threads 1 -o /gs7/users/user/eggnog-mapper/test.emapper.hits  --sensitive --iterate -e 0.001 --top 3  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp
Functional annotation of hits...
0 3.5762786865234375e-06 0.00 q/s (% mem usage: 3.90, % mem avail: 96.09)
2 2.5320394039154053 0.79 q/s (% mem usage: 3.90, % mem avail: 96.08)
Done
   /data/user/eggnog-mapper/test.emapper.hits
   /data/user/eggnog-mapper/test.emapper.seed_orthologs
   /data/user/eggnog-mapper/test.emapper.annotations
...

```

Alternatively, you can create a diamond database on your own:

```

[user@cn3200 ~]$ **mkdir data** 
[user@cn3200 ~]$ **export EGGNOG\_DATA\_DIR=./data** 
[user@cn3200 ~]$ **create\_dbs.py -m diamond --dbname bacteria --taxa Bacteria**

```

This will create a bacteria.dmnd diamond database in the directory specified in EGGNOG\_DATA\_DIR environment variable.





