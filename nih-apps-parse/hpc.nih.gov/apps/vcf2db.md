

document.querySelector('title').textContent = 'Vcf2db on Biowulf';
Vcf2db on Biowulf


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



Previously (and currently), gemini kept a bunch of vetted annotations along with the gemini install and annotated an incoming VCF with those annotations as it was loaded into gemini. This is nice for users but by de-coupling the annotation from the loading, we have more flexiblility.



Documentation
* [vcf2db Main Site](https://github.com/quinlan-lab/vcf2db)


Important Notes
* Module Name: vcf2db (see [the modules page](/apps/modules.html) for more information)
* singlethreaded app
* Example files in /usr/local/apps/vcf2db/tests



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

[user@cn3144 ~]$ **cp /usr/local/apps/vcf2db/tests/test.\* .**

[user@cn3144 ~]$ **module load vcf2db**

[user@cn3144 ~]$ **vcf2db.py -h**
usage: 
Take a VCF and create a gemini compatible database

       [-h] [--a-ok A_OK] [-e INFO_EXCLUDE] [--impacts-field IMPACTS_FIELD]
       [--legacy-compression]
       [--expand {gt_ref_depths,gt_types,gt_depths,gt_alt_depths,gt_quals,gt_alt_freqs}]
       VCF ped db

positional arguments:
  VCF
  ped
  db

optional arguments:
  -h, --help            show this help message and exit
  --a-ok A_OK           list of info names to include even with Number=A (will
                        error if they have > 1 value
  -e INFO_EXCLUDE, --info-exclude INFO_EXCLUDE
                        don't save this field to the database. May be
                        specified multiple times.
  --impacts-field IMPACTS_FIELD
                        this field should be propagated to the variant_impacts
                        table. by default, only CSQ/EFF/ANN fields are added.
                        the field can be suffixed with a type of ':i' or ':f'
                        to indicate int or float to override the default of
                        string. e.g. AF:f
  --legacy-compression
  --expand {gt_ref_depths,gt_types,gt_depths,gt_alt_depths,gt_quals,gt_alt_freqs}
                        sample columns to expand into their own tables


[user@cn3144 ~]$ **vcf2db.py test.vcf test.ped test.db**
skipping 'AC' because it has Number=A
skipping 'AF' because it has Number=A
pedigree notice: '1_dad' is dad but has unknown sex. Setting to male
pedigree notice: '1_mom' is mom but has unknown sex. Setting to female
pedigree notice: '2_dad' is dad but has unknown sex. Setting to male
pedigree notice: '2_mom' is mom but has unknown sex. Setting to female
pedigree notice: '3_dad' is dad but has unknown sex. Setting to male
pedigree notice: '3_mom' is mom but has unknown sex. Setting to female
/usr/local/Anaconda/envs_app/vcf2db/2017.12.11/lib/python2.7/site-packages/sqlalchemy/sql/sqltypes.py:226: SAWarning: Unicode type received non-unicode bind param value '1'. (this warning may be suppressed after 10 occurrences)
  (util.ellipses_string(value),))
/usr/local/Anaconda/envs_app/vcf2db/2017.12.11/lib/python2.7/site-packages/sqlalchemy/sql/sqltypes.py:226: SAWarning: Unicode type received non-unicode bind param value 'aaa'. (this warning may be suppressed after 10 occurrences)
  (util.ellipses_string(value),))
/usr/local/Anaconda/envs_app/vcf2db/2017.12.11/lib/python2.7/site-packages/sqlalchemy/sql/sqltypes.py:226: SAWarning: Unicode type received non-unicode bind param value 'bbb'. (this warning may be suppressed after 10 occurrences)
  (util.ellipses_string(value),))
/usr/local/Anaconda/envs_app/vcf2db/2017.12.11/lib/python2.7/site-packages/sqlalchemy/sql/sqltypes.py:226: SAWarning: Unicode type received non-unicode bind param value '2'. (this warning may be suppressed after 10 occurrences)
  (util.ellipses_string(value),))
/usr/local/Anaconda/envs_app/vcf2db/2017.12.11/lib/python2.7/site-packages/sqlalchemy/sql/sqltypes.py:226: SAWarning: Unicode type received non-unicode bind param value 'ccc'. (this warning may be suppressed after 10 occurrences)
  (util.ellipses_string(value),))
/usr/local/Anaconda/envs_app/vcf2db/2017.12.11/lib/python2.7/site-packages/sqlalchemy/sql/sqltypes.py:226: SAWarning: Unicode type received non-unicode bind param value 'ddd'. (this warning may be suppressed after 10 occurrences)
  (util.ellipses_string(value),))
/usr/local/Anaconda/envs_app/vcf2db/2017.12.11/lib/python2.7/site-packages/sqlalchemy/sql/sqltypes.py:226: SAWarning: Unicode type received non-unicode bind param value 'eee'. (this warning may be suppressed after 10 occurrences)
  (util.ellipses_string(value),))
/usr/local/Anaconda/envs_app/vcf2db/2017.12.11/lib/python2.7/site-packages/sqlalchemy/sql/sqltypes.py:226: SAWarning: Unicode type received non-unicode bind param value 'fff'. (this warning may be suppressed after 10 occurrences)
  (util.ellipses_string(value),))
6 variant_impacts:23  effects time: 0.0 chunk time:0.1  43.33 variants/second
indexing ... finished in 0.0 seconds...
total time: in 0.2 seconds...

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. vcf2db.sh). For example:



```

#!/bin/bash
module load vcf2db
vcf2db.py test.vcf test.ped test.db

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch vcf2db.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. vcf2db.swarm). For example:



```

vcf2db.py test1.vcf test1.ped test1.db
vcf2db.py test2.vcf test2.ped test2.db
vcf2db.py test3.vcf test3.ped test3.db

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f vcf2db.swarm [-g 4] [-t 2] --module vcf2db
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module vcf2db Loads the vcf2db module for each subjob in the swarm 
 | |
 | |
 | |








