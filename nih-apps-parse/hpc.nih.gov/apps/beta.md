

document.querySelector('title').textContent = 'BETA: Binding and Expression Target Analysis ';
**BETA: Binding and Expression Target Analysis** 


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



Binding and expression target analysis (BETA) is a software package 
that integrates ChIP-seq of TFs or chromatin regulators 
with differential gene expression data to infer direct target genes. 
The combination of ChIP-seq and transcriptome analysis is a compelling approach to unravel the regulation of gene expression. 



### Reference:


* Su Wang, Hanfei Sun, Jian Ma, Chongzhi Zang, Chenfei Wang, Juan Wang, Qianzi Tang, Clifford A Meyer, Yong Zhang and X Shirley Liu   

*Target analysis by integration of transcriptome and ChIP-seq data with BETA*   

[Nature Protocols](https://www.nature.com/articles/nprot.2013.150) 2013, **8**(12), p.2502-2515.


Documentation
* [Citrome BETA Home page](http://cistrome.org/BETA)
* [BETA Github\_page](https://github.com/suwangbio/BETA)


Important Notes
* Module Name: BETA (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **BETA\_HOME**  BETA installation directory
	+ **BETA\_BIN**       BETA executable directory
	+ **BETA\_DOC**       BETA documentation directory
	+ **BETA\_DATA**    BETA test data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive** 
[user@cn3111 ~]$ **module load BETA**
[user@cn3111 ~]$ **BETA -h**
usage: BETA [-h] [--version] {basic,plus,minus} ...

BETA --- Binding Expression Target Analysis

positional arguments:
  {basic,plus,minus}  sub-command help
    basic             Main BETA Function: Transcription factors direct targets
                      prediction.
    plus              The super beta beta can not only do TF's direct targets
                      and active or repressive prediction, but also do the
                      motif analysis with the target regions
    minus             Find Target Genes with only binding data: regulatiry
                      potential score

optional arguments:
  -h, --help          show this help message and exit
  --version           show program's version number and exit

For command line options of each command, type: BETA COMMAND -h

```

Download sample data:

```

[user@cn3111 ~]$ **cp $BETA\_DATA/\* .** 

```

Link a human genome sequence in FASTA format to the current folder:

```

[user@cn3111 ~]$ **ln -s /fdb/igenomes/Homo\_sapiens/UCSC/hg19/hg19.fa**

```

BETA Basic will do the factor function prediction and direct target detecting:

```

[user@cn3111 ~]$ **BETA basic -p 3656\_peaks.bed -e AR\_diff\_expr.xls -k LIM -g hg19 --da 500 -n basic --info 1,2,6** 
[09:09:02] Argument List: 
[09:09:02] Name = basic
[09:09:02] Peak File = 3656_peaks.bed
[09:09:02] Top Peaks Number = 10000
[09:09:02] Distance = 100000 bp
[09:09:02] Genome = hg19
[09:09:02] Expression File = AR_diff_expr.xls
[09:09:02] Expression Type = MicroArray, LIMMA result
[09:09:02] Number of differential expressed genes = 500.0
[09:09:02] Differential expressed gene FDR Threshold = 1
[09:09:02] Up/Down Prediction Cutoff = 0.001000
[09:09:02] Function prediction based on regulatory potential
[09:09:02] Check 3656_peaks.bed successfully!
[09:09:02] #ID	logFC	AveExpr	t	P.Value	adj.P.Val	B
 is the header of the expression file
[09:09:02] Checking the differential expression infomation...
[09:09:02] Take the first line with Differential Information as an example: NR_045762_at3.16711734	9.140369116	35.91057535	6.99E-11	4.18E-07	14.13456018

[09:09:02] Differential Expression file format successful passed
[09:09:02] You do not like filter peak by CFCT boundary, it will be filtered only by the distance
[09:09:02] Read file <3656_peaks.bed> OK! All <7059> peaks.
[09:09:09] Process <50801> genes
[09:09:16] Finished! Preliminary results saved into temporary file: 
[15841, 18956]
[09:09:16] Genes were seprated to two parts: up regulated and down regulated.
cut: write error: Broken pipe
cut: write error: Broken pipe
[09:09:17] Prepare file for the Up/Down Test
null device 
 1 
[09:09:20] Finished, Find the result in basic\_function\_prediction.pdf
[09:09:20] Get the Rank Product of the "upregulate" genes
['"upregulate"']
[09:10:00] pick out the peaks 100000bp around the selected genes
...
total time: 0:0:58 


```

BETA Plus will do TF active and repressive function prediction, direct targets
detecting and motif analysis in target regions:

```

[user@cn3111 ~]$ **BETA plus -p 3656\_peaks.bed -e AR\_diff\_expr.xls -k LIM -g hg19 --gs hg19.fa** 
[08:38:05] Argument List: 
[08:38:05] Name = NA
[08:38:05] Peak File = 3656_peaks.bed
[08:38:05] Top Peaks Number = 10000
[08:38:05] Distance = 100000 bp
[08:38:05] Genome = hg19
[08:38:05] Expression File = AR_diff_expr.xls
[08:38:05] Genome Sequence fasta formated data = hg19.fa
[08:38:05] Expression Type = MicroArray, LIMMA result
[08:38:05] Number of differential expressed genes = 0.5
[08:38:05] Differential expressed gene FDR Threshold = 1
[08:38:05] Up/Down Prediction Cutoff = 0.001000
[08:38:05] Function prediction based on regulatory potential
[08:38:05] Check 3656_peaks.bed successfully!
[08:38:05] #ID	logFC	AveExpr	t	P.Value	adj.P.Val	B
 is the header of the expression file
[08:38:05] Checking the differential expression infomation...
[08:38:05] Take the first line with Differential Information as an example: NR_045762_at3.16711734	9.140369116	35.91057535	6.99E-11	4.18E-07	14.13456018

[08:38:05] Differential Expression file format successful passed
[08:38:05] You do not like filter peak by CFCT boundary, it will be filtered only by the distance
[08:38:06] Read file <3656_peaks.bed> OK! All <7059> peaks.
[08:38:13] Process <50801> genes
[08:38:20] Finished! Preliminary results saved into temporary file: 
[15841, 18956]
[08:38:20] Genes were seprated to two parts: up regulated and down regulated.
cut: write error: Broken pipe
cut: write error: Broken pipe
[08:38:20] Prepare file for the Up/Down Test
null device 
 1 
[08:39:05] Finished, Find the result in NA\_function\_prediction.pdf
[08:39:05] Get the Rank Product of the "upregulate" genes
['"upregulate"']
[08:39:35] pick out the peaks 100000bp around the selected genes
...
[08:39:36] get three regions of every peak around the up genes
['up\_middle.bed', 'up\_left.bed', 'up\_right.bed']
 [08:40:42] get the fasta format sequence data of the three regions of up
[08:40:42] get three regions of every peak around the non genes
['non\_middle.bed', 'non\_left.bed', 'non\_right.bed']
[08:41:46] get the fasta format sequence data of the three regions of non
[08:41:46] run mis to do the known motif scan with cistrome motif database
Calculating Background..
A or T: 161685
C or G: 131715
Others: 0
 Calculating Background..
A or T: 161503
C or G: 130430
Others: 0
Calculating Background..
A or T: 159938
C or G: 131995
Others: 0
Calculating Background..
A or T: 236519
C or G: 149281
Others: 0
Calculating Background..
A or T: 236826
C or G: 147045
Others: 0
Calculating Background..
A or T: 235215
C or G: 148656
Others: 0
[08:43:51] T statistic test performed here to do the significance testing of every motif
[09:05:12] motif result will be in html format
...
Success parser from table.
[09:05:20] Done: find motif result in beatmotif.html file
...
total time: 0:27:14 


```

End the interactive session:

```

[user@cn3111 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





