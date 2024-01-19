

document.querySelector('title').textContent = 'Biom-format on Biowulf';
Biom-format on Biowulf


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



The biom-format python package contains a command line tool to manipulate
and convert [biom format](http://biom-format.org/) files. It
also includes an API for programatically manipulating biom files. It is
therefore installed as both in independent application and as part of the
python environments.



### References:


* Daniel McDonald *et al.*. *The Biological Observation
 Matrix (BIOM) format or: how I learned to stop worrying and love
 the ome-ome*. GigaScience 1 (2012), doi:10.1186/2047-217X-1-7.
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/23587224) | 
 [PMC](http://www.ncbi.nlm.nih.gov/pmc/articles/pmid/23587224/) | 
 [Journal](http://www.gigasciencejournal.com/content/1/1/7)


Documentation
* [biom-format Main Site](http://biom-format.org/)


Important Notes
* Module Name: biom-format (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded app
* Example files in /usr/local/apps/biom-format/TEST\_DATA



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session (user input in **bold**:



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load biom-format**

[user@cn3144 ~]$ **TD=/usr/local/apps/biom-format/TEST\_DATA**

[user@cn3144 ~]$ **biom head -i $TD/phinch\_testdata.biom**
# Constructed from biom file
#OTU ID 0.IntakeWater.1 0.IntakeWater.3 0.IntakeWater.2 27.WaterCoralpond.2 
228057  3.0     6.0     3.0     1.0     0.0
988537  0.0     0.0     0.0     2.0     0.0
89370   0.0     0.0     0.0     1.0     1.0
2562097 0.0     0.0     0.0     0.0     0.0
256904  8.0     0.0     1.0     32.0    58.0

[user@cn3144 ~]$ **biom convert -i $TD/phinch\_testdata.biom -o test.biom --to-hdf5**

[user@cn3144 ~]$ **module load hdf5**

[user@cn3144 ~]$ **h5ls test.biom**
observation              Group
sample                   Group

[user@cn3144 ~]$ **h5ls test.biom/sample**
group-metadata           Group
ids                      Dataset {95}
matrix                   Group
metadata                 Group

[user@cn3144 ~]$ **biom summarize-table -i test.biom**
Num samples: 95
Num observations: 67900
Total count: 10223009
Table density (fraction of non-zero values): 0.056

Counts/sample summary:
 Min: 16.0
 Max: 1106184.0
 Median: 65463.000
 Mean: 107610.621
 Std. dev.: 164313.408
 Sample Metadata Categories: description; alkalinity; material; ammonium; nitrite; 
  LinkerPrimerSequence; sulfide; BarcodeName; InternalCode; temp; 
  collection_date; BarcodeSequence; salinity; phosphate; ReverseBarcode; nitrate; 
  ph; ReverseName; ReversePrimerSequence; Hardness; diss_oxygen
 Observation Metadata Categories: taxonomy

Counts/sample detail:
0.WipesKoipondLgWaterfall.1: 16.0
0.WipesKoipondLFilter.1: 18.0
[...snip...]

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

When loading one of the python 2.7 modules, the `biom`
command line tool will also become available (though the version may vary
over time), as will the API. For example




```

[user@cn3144 ~]$ **module load python**

[user@cn3144 ~]$ **which biom**
/usr/local/Anaconda/envs/py2.7.9/bin/biom

[user@cn3144 ~]$ **python**
Python 2.7.9 |Continuum Analytics, Inc.| (default, Apr 14 2015, 12:54:25)
[GCC 4.4.7 20120313 (Red Hat 4.4.7-1)] on linux2
Type "help", "copyright", "credits" or "license" for more information.
Anaconda is brought to you by Continuum Analytics.
Please check out: http://continuum.io/thanks and https://binstar.org

>>> **import biom**

>>> **table = biom.load\_table("test.biom")**

>>> **table**
67900 x 95 <class 'biom.table.Table'> with 360783 nonzero entries (5% dense)

```








