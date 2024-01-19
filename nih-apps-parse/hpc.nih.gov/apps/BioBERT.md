

document.querySelector('title').textContent = 'BioBERT: a biomedical language representation model designed for biomedical text mining tasks';
**BioBERT: a biomedical language representation model   
designed for biomedical text mining tasks**


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



BioBERT is a biomedical language representation model 
designed for biomedical text mining tasks 
such as biomedical named entity recognition, relation extraction, question answering, etc. 



### References:


* Jinhyuk Lee, Wonjin Yoon, Sungdong Kim, Donghyeon Kim, Sunkyu Kim, Chan Ho So and Jaewoo Kang,   

*BioBERT: a pre-trained biomedical language representation model for biomedical text mining*   

[Bioinformatics](https://academic.oup.com/bioinformatics/article/36/4/1234/5566506) (2020),
**36**(4), 1234–1240. doi: 10.1093/bioinformatics/btz682.


Documentation
* [BioBERT Github page](https://github.com/dmis-lab/biobert)
* [BioBERT Tutorial](https://towardsdatascience.com/tagging-genes-and-proteins-with-biobert-c7b04fc6eb4f)


Important Notes
* Module Name: BioBERT (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **BIOBERT\_HOME**  installation directory
	+ **BIOBERT\_BIN**  executable directory
	+ **BIOBERT\_DIR**  source code directory
	+ **BIOBERT\_DATA**  configuration files and models directory
	+ **BIOBERT\_BENCHMARKS**  benchmark datasets directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=16g --gres=gpu:p100:1,lscratch:10 -c4** 
[user@cn3104 ~]$ **module load biobert**
[+] Loading singularity  3.10.5  on cn3104
[+] Loading BioBERT v20200409

```

Copy sample data to your current folder: and run one of BioBert scripts:

```

[user@cn3104 ~]$ **cp $BIOBERT\_DATA/\* .**
[user@cn3104 ~]$ **export BIOBERT\_DIR=.** 
[user@cn3104 ~]$ **cp -r $BIOBERT\_BENCHMARKS/NCBI-disease .**
[user@cn3104 ~]$ **export NER\_DIR=./NCBI-disease** 
[user@cn3104 ~]$ **export OUTPUT\_DIR=./ner\_outputs** 
[user@cn3104 ~]$ **mkdir -p $OUTPUT\_DIR** 

[user@cn3104 ~]$ **python $BIOBERT\_SRC/run\_ner.py --do\_train=true --do\_eval=true --vocab\_file=$BIOBERT\_DIR/vocab.txt --bert\_config\_file=$BIOBERT\_DIR/bert\_config.json --init\_checkpoint=$BIOBERT\_DIR/model.ckpt-1000000 --num\_train\_epochs=10.0 --data\_dir=$NER\_DIR --output\_dir=$OUTPUT\_DIR**
...
INFO:tensorflow:Writing example 0 of 8230
I1012 11:53:03.600488 23456244467520 run_ner.py:303] Writing example 0 of 8230
INFO:tensorflow:*** Example ***
I1012 11:53:03.601105 23456244467520 run_ner.py:276] *** Example ***
INFO:tensorflow:guid: train-0
I1012 11:53:03.601165 23456244467520 run_ner.py:277] guid: train-0
INFO:tensorflow:tokens: I ##dent ##ification of AP ##C ##2 , a ho ##mo ##logue of the ad ##eno ##mat ##ous p ##oly ##po ##sis co ##li t ##umour suppress ##or .
I1012 11:53:03.601217 23456244467520 run_ner.py:279] tokens: I ##dent ##ification of AP ##C ##2 , a ho ##mo ##logue of the ad ##eno ##mat ##ous p ##oly ##po ##sis co ##li t ##umour suppress ##or .
INFO:tensorflow:input_ids: 101 146 11951 5783 1104 10997 1658 1477 117 170 16358 3702 12733 1104 1103 8050 26601 21943 2285 185 23415 5674 4863 1884 2646 189 27226 17203 1766 119 102 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
I1012 11:53:03.601279 23456244467520 run_ner.py:280] input_ids: 101 146 11951 5783 1104 10997 1658 1477 117 170 16358 3702 12733 1104 1103 8050 26601 21943 2285 185 23415 5674 4863 1884 2646 189 27226 17203 1766 119 102 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
INFO:tensorflow:input_mask: 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
I1012 11:53:03.601342 23456244467520 run_ner.py:281] input_mask: 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
INFO:tensorflow:segment_ids: 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
I1012 11:53:03.601401 23456244467520 run_ner.py:282] segment_ids: 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
INFO:tensorflow:label_ids: 5 3 4 4 3 3 4 4 3 3 3 4 4 3 3 1 4 4 4 2 4 4 4 2 4 2 4 3 4 3 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
I1012 11:53:03.601460 23456244467520 run_ner.py:283] label_ids: 5 3 4 4 3 3 4 4 3 3 3 4 4 3 3 1 4 4 4 2 4 4 4 2 4 2 4 3 4 3 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
INFO:tensorflow:*** Example ***
I1012 11:53:03.602299 23456244467520 run_ner.py:276] *** Example ***
INFO:tensorflow:guid: train-1
I1012 11:53:03.602358 23456244467520 run_ner.py:277] guid: train-1
INFO:tensorflow:tokens: The ad ##eno ##mat ##ous p ##oly ##po ##sis co ##li ( AP ##C ) t ##umour - suppress ##or protein controls the W ##nt signalling pathway by forming a complex with g ##ly ##co ##gen s ##ynth ##ase kinase 3 ##bet ##a ( G ##S ##K - 3 ##bet ##a )
I1012 11:53:03.602412 23456244467520 run_ner.py:279] tokens: The ad ##eno ##mat ##ous p ##oly ##po ##sis co ##li ( AP ##C ) t ##umour - suppress ##or protein controls the W ##nt signalling pathway by forming a complex with g ##ly ##co ##gen s ##ynth ##ase kinase 3 ##bet ##a ( G ##S ##K - 3 ##bet ##a )
INFO:tensorflow:input_ids: 101 1109 8050 26601 21943 2285 185 23415 5674 4863 1884 2646 113 10997 1658 114 189 27226 118 17203 1766 4592 7451 1103 160 2227 25498 13548 1118 5071 170 2703 1114 176 1193 2528 4915 188 26588 6530 24779 124 16632 1161 113 144 1708 2428 118 124 16632 1161 114 102 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
I1012 11:53:03.602473 23456244467520 run_ner.py:280] input_ids: 101 1109 8050 26601 21943 2285 185 23415 5674 4863 1884 2646 113 10997 1658 114 189 27226 118 17203 1766 4592 7451 1103 160 2227 25498 13548 1118 5071 170 2703 1114 176 1193 2528 4915 188 26588 6530 24779 124 16632 1161 113 144 1708 2428 118 124 16632 1161 114 102 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
INFO:tensorflow:input_mask: 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
I1012 11:53:03.602534 23456244467520 run_ner.py:281] input_mask: 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
INFO:tensorflow:segment_ids: 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
I1012 11:53:03.602594 23456244467520 run_ner.py:282] segment_ids: 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
INFO:tensorflow:label_ids: 5 3 1 4 4 4 2 4 4 4 2 4 2 2 4 2 2 4 3 3 4 3 3 3 3 4 3 3 3 3 3 3 3 3 4 4 4 3 4 4 3 3 4 4 3 3 4 4 3 3 4 4 3 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
I1012 11:53:03.602653 23456244467520 run_ner.py:283] label_ids: 5 3 1 4 4 4 2 4 4 4 2 4 2 2 4 2 2 4 3 3 4 3 3 3 3 4 3 3 3 3 3 3 3 3 4 4 4 3 4 4 3 3 4 4 3 3 4 4 3 3 4 4 3 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
INFO:tensorflow:*** Example ***
I1012 11:53:03.603083 23456244467520 run_ner.py:276] *** Example ***
INFO:tensorflow:guid: train-2
I1012 11:53:03.603142 23456244467520 run_ner.py:277] guid: train-2
INFO:tensorflow:tokens: , a ##xin / conduct ##in and beta ##cate ##nin .
I1012 11:53:03.603189 23456244467520 run_ner.py:279] tokens: , a ##xin / conduct ##in and beta ##cate ##nin .
INFO:tensorflow:input_ids: 101 117 170 16594 120 5880 1394 1105 11933 20127 10430 119 102 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
I1012 11:53:03.603248 23456244467520 run_ner.py:280] input_ids: 101 117 170 16594 120 5880 1394 1105 11933 20127 10430 119 102 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
INFO:tensorflow:input_mask: 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
I1012 11:53:03.603308 23456244467520 run_ner.py:281] input_mask: 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
INFO:tensorflow:segment_ids: 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
I1012 11:53:03.603368 23456244467520 run_ner.py:282] segment_ids: 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
INFO:tensorflow:label_ids: 5 3 3 4 3 3 4 3 3 4 4 3 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
I1012 11:53:03.603425 23456244467520 run_ner.py:283] label_ids: 5 3 3 4 3 3 4 3 3 4 4 3 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
INFO:tensorflow:*** Example ***
I1012 11:53:03.603892 23456244467520 run_ner.py:276] *** Example ***
INFO:tensorflow:guid: train-3
I1012 11:53:03.603949 23456244467520 run_ner.py:277] guid: train-3
INFO:tensorflow:tokens: Complex formation induce ##s the rapid degradation of beta ##cate ##nin .
I1012 11:53:03.603994 23456244467520 run_ner.py:279] tokens: Complex formation induce ##s the rapid degradation of beta ##cate ##nin .
INFO:tensorflow:input_ids: 101 9974 3855 21497 1116 1103 6099 18126 1104 11933 20127 10430 119 102 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
I1012 11:53:03.604052 23456244467520 run_ner.py:280] input_ids: 101 9974 3855 21497 1116 1103 6099 18126 1104 11933 20127 10430 119 102 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
INFO:tensorflow:input_mask: 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
I1012 11:53:03.604111 23456244467520 run_ner.py:281] input_mask: 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
INFO:tensorflow:segment_ids: 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
I1012 11:53:03.604169 23456244467520 run_ner.py:282] segment_ids: 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
INFO:tensorflow:label_ids: 5 3 3 3 4 3 3 3 3 3 4 4 3 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
I1012 11:53:03.604226 23456244467520 run_ner.py:283] label_ids: 5 3 3 3 4 3 3 3 3 3 4 4 3 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
INFO:tensorflow:*** Example ***
I1012 11:53:03.604921 23456244467520 run_ner.py:276] *** Example ***
INFO:tensorflow:guid: train-4
I1012 11:53:03.604977 23456244467520 run_ner.py:277] guid: train-4
INFO:tensorflow:tokens: In co ##lon car ##cin ##oma cells , loss of AP ##C leads to the accumulation of beta ##cate ##nin in the nucleus , where it binds to and activate ##s the T ##c ##f - 4 transcription factor
I1012 11:53:03.605027 23456244467520 run_ner.py:279] tokens: In co ##lon car ##cin ##oma cells , loss of AP ##C leads to the accumulation of beta ##cate ##nin in the nucleus , where it binds to and activate ##s the T ##c ##f - 4 transcription factor
INFO:tensorflow:input_ids: 101 1130 1884 4934 1610 16430 7903 3652 117 2445 1104 10997 1658 4501 1106 1103 23168 1104 11933 20127 10430 1107 1103 14297 117 1187 1122 23126 1106 1105 23162 1116 1103 157 1665 2087 118 125 15416 5318 102 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
I1012 11:53:03.605087 23456244467520 run_ner.py:280] input_ids: 101 1130 1884 4934 1610 16430 7903 3652 117 2445 1104 10997 1658 4501 1106 1103 23168 1104 11933 20127 10430 1107 1103 14297 117 1187 1122 23126 1106 1105 23162 1116 1103 157 1665 2087 118 125 15416 5318 102 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
INFO:tensorflow:input_mask: 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
I1012 11:53:03.605145 23456244467520 run_ner.py:281] input_mask: 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
INFO:tensorflow:segment_ids: 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
I1012 11:53:03.605203 23456244467520 run_ner.py:282] segment_ids: 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
INFO:tensorflow:label_ids: 5 3 1 4 2 4 4 3 3 3 3 3 4 3 3 3 3 3 3 4 4 3 3 3 3 3 3 3 3 3 3 4 3 3 4 4 3 3 3 3 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
I1012 11:53:03.605262 23456244467520 run_ner.py:283] label_ids: 5 3 1 4 2 4 4 3 3 3 3 3 4 3 3 3 3 3 3 4 4 3 3 3 3 3 3 3 3 3 3 4 3 3 4 4 3 3 3 3 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
INFO:tensorflow:Writing example 5000 of 8230
...

```

Run BioBERT on the QA dataset:

```

[user@cn3104 ~]$ **cp $BIOBERT\_DATA/\* .**
[user@cn3104 ~]$ **export BIOBERT\_DIR=.** 
[user@cn3104 ~]$ **cp -r $BIOBERT\_BENCHMARKS/BioASQ .**
[user@cn3104 ~]$ **export QA\_DIR=BioASQ**
[user@cn3104 ~]$ **export OUTPUT\_DIR=./qa\_outputs**
[user@cn3104 ~]$ **mkdir -p $OUTPUT\_DIR**
[user@cn3104 ~]$ **python $BIOBERT\_SRC/run\_qa.py --do\_train=True --do\_predict=True --vocab\_file=$BIOBERT\_DIR/vocab.txt --bert\_config\_file=$BIOBERT\_DIR/bert\_config.json --init\_checkpoint=$BIOBERT\_DIR/model.ckpt-1000000 --max\_seq\_length=384 --train\_batch\_size=12 --learning\_rate=5e-6 --doc\_stride=128 --num\_train\_epochs=5.0 --do\_lower\_case=False --train\_file=$QA\_DIR/BioASQ-train-factoid-4b.json --predict\_file=$QA\_DIR/BioASQ-test-factoid-4b-1.json --output\_dir=$OUTPUT\_DIR** 
...

Exit the application:   
      

```

[user@cn3104 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





```


