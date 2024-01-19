

document.querySelector('title').textContent = 'Bartender: fast and accurate counting of barcode reads';
**Bartender: fast and accurate counting of barcode reads**


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



Bartender is an accurate clustering algorithm to detect barcodes
and their abundances from raw next-generation sequencing data.
In contrast with existing methods that cluster based on sequence similarity alone,
Bartender uses a modified two-sample proportion test that also considers cluster size. This modification
results in higher accuracy and lower rates of under- and over-clustering artifacts.



### References:


* Lu Zhao, Zhimin Liu, Sasha F. Levy and Song Wu  

*Bartender: a fast and accurate clustering algorithm to count barcode reads*    

[Bioinformatics](https://academic.oup.com/bioinformatics/article/34/5/739/4562326?login=true) **34**(5),
739–747 (2018).


Documentation
* [Bartender GitHub page](https://github.com/LaoZZZZZ/bartender-1.1)


Important Notes
* Module Name: Bartender (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **BARTENDER\_HOME**  installation directory
	+ **BARTENDER\_BIN**       executable directory
	+ **BARTENDER\_SRC**       source code directory
	+ **BARTENDER\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive** 
[user@cn3200 ~]$**module load bartender** 
[+] Loading singularity  3.8.2  
[+] Loading bartender  1.1
[user@cn3200 ~]$ **cp $BARTENDER\_DATA/\* .**
[user@cn3200 ~]$ **bartender\_extractor\_com -f 2M\_test.fq -o 2M\_extracted -q "?" -p "TACC[4-7]AA[4-7]AA[4-7]TT[4-7]ATAA" -m 2**
    Running bartender extractor
    bartender_extractor 2M_test.fq 2M_extracted 63 "(TAC.|TA.C|T.CC|.ACC)([ATCGN]{4,7})(AA)([ATCGN]{4,7})(AA)([ATCGN]{4,7})(TT)([ATCGN]{4,7})(ATA.|AT.A|A.AA|.TAA)" TACC ATAA 9 1
    Totally there are 1000 reads in 2M_test.fq file!
    Totally there are 976 valid barcodes from 2M_test.fq file
    Totally there are 924 valid barcodes whose quality pass the quality condition
    The estimated sequence error from the prefix and suffix parts is 0.00153689

[user@cn3200 ~]$ **bartender\_single\_com -f 2M\_extracted\_barcode.txt -o 2M\_barcode -d 3**
Running bartender
Loading barcodes from the file
It takes 00:00:00 to load the barcodes from 2M_extracted_barcode.txt
Shortest barcode length: 18
Longest barcode length: 26
Start to group barcode with length 18
Using two sample unpooled test
Transforming the barcodes into seed clusters
Initial number of unique reads:  1
The distance threshold is 3
Identified 1 barcodes with length 18
Start to group barcode with length 19
Using two sample unpooled test
Transforming the barcodes into seed clusters
Initial number of unique reads:  7
The distance threshold is 3
Clustering iteration 1
Clustering iteration 2
Clustering iteration 3
Clustering iteration 4
Identified 7 barcodes with length 19
Start to group barcode with length 20
Using two sample unpooled test
Transforming the barcodes into seed clusters
Initial number of unique reads:  221
The distance threshold is 3
Clustering iteration 1
Clustering iteration 2
Clustering iteration 3
Clustering iteration 4
Identified 221 barcodes with length 20
Start to group barcode with length 21
Using two sample unpooled test
Transforming the barcodes into seed clusters
Initial number of unique reads:  0
The distance threshold is 3
Identified 0 barcodes with length 21
Start to group barcode with length 22
Using two sample unpooled test
Transforming the barcodes into seed clusters
Initial number of unique reads:  0
The distance threshold is 3
Identified 0 barcodes with length 22
Start to group barcode with length 23
Using two sample unpooled test
Transforming the barcodes into seed clusters
Initial number of unique reads:  0
The distance threshold is 3
Identified 0 barcodes with length 23
Start to group barcode with length 24
Using two sample unpooled test
Transforming the barcodes into seed clusters
Initial number of unique reads:  0
The distance threshold is 3
Identified 0 barcodes with length 24
Start to group barcode with length 25
Using two sample unpooled test
Transforming the barcodes into seed clusters
Initial number of unique reads:  0
The distance threshold is 3
Identified 0 barcodes with length 25
Start to group barcode with length 26
Using two sample unpooled test
Transforming the barcodes into seed clusters
Initial number of unique reads:  2
The distance threshold is 3
Clustering iteration 1
Clustering iteration 2
Clustering iteration 3
Identified 2 barcodes with length 26
The clustering process takes 00:00:00
Start to dump clusters to file with prefix 2M_barcode
Start to remove pcr effects
***(Overall error rate estimated from the clustering result)***
Total number of clusters after removing PCR effects: 231
Could not find any high quality clusters(max entropy < 0.33, cluster size > 20 ) to estimate the error rate!
Please use the result cautiously(Better to check the input data)!

The estimated error rate is 0
The overall running time 00:00:00 seconds.

```


End the interactive session:

```

[user@cn3200 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





