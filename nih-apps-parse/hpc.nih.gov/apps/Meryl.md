

document.querySelector('title').textContent = 'Meryl: a genomic k-mer counter (and sequence utility) with nice features';
**Meryl: a genomic k-mer counter (and sequence utility) with nice features**


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



Meryl is the k-mer counter. It is built into the Celera assembler 
and is also available as a stand-alone application. 
Meryl uses a sorting-based approach that sorts the k-mers in lexicographical order. 



### References:


* Jason R. Miller, Arthur L. Delcher, Sergey Koren, Eli Venter, Brian P. Walenz, Anushka Brownley, Justin Johnson, Kelvin Li, Clark Mobarry, Granger Sutton.
*Aggressive assembly of pyrosequencing reads with mates,*   

[Bioinformatics](https://academic.oup.com/bioinformatics/article/24/24/2818/197033) 2008, vol. 24 (pg. 2818-2824).


Documentation
* [Meryl GitHub page](https://github.com/marbl/meryl)


Important Notes
* Module Name: Meryl (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Implemented as a Singularity container
* Unusual environment variables set
	+ **MERYL\_HOME**  installation directory
	+ **MERYL\_HOME**  executable directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=4g**
[user@cn3316 ~]$ **module load meryl**
[+] Loading Meryl 0.0  ...
[user@cn3316 ~]$ **meryl -h**
meryl -h
usage: meryl ...

  A meryl command line is formed as a series of commands and files, possibly
  grouped using square brackets.  Each command operates on the file(s) that
  are listed after it.

  COMMANDS:

    print                display kmers on the screen as 'kmercount'. accepts exactly one input.

 count Count the occurrences of canonical kmers in the input. must have 'output' specified.
 count-forward Count the occurrences of forward kmers in the input. must have 'output' specified.
 count-reverse Count the occurrences of reverse kmers in the input. must have 'output' specified.
 k= create mers of size K bases (mandatory).
 n= expect N mers in the input (optional; for precise memory sizing).
 memory=M use no more than (about) M GB memory.
 threads=T use no more than T threads.

 less-than N return kmers that occur fewer than N times in the input. accepts exactly one input.
 greater-than N return kmers that occur more than N times in the input. accepts exactly one input.
 equal-to N return kmers that occur exactly N times in the input. accepts exactly one input.
 not-equal-to N return kmers that do not occur exactly N times in the input. accepts exactly one input.

 increase X add X to the count of each kmer.
 decrease X subtract X from the count of each kmer.
 multiply X multiply the count of each kmer by X.
 divide X divide the count of each kmer by X.
 modulo X set the count of each kmer to the remainder of the count divided by X.

 union return kmers that occur in any input, set the count to the number of inputs with this kmer.
 union-min return kmers that occur in any input, set the count to the minimum count
 union-max return kmers that occur in any input, set the count to the maximum count
 union-sum return kmers that occur in any input, set the count to the sum of the counts

 intersect return kmers that occur in all inputs, set the count to the count in the first input.
 intersect-min return kmers that occur in all inputs, set the count to the minimum count.
 intersect-max return kmers that occur in all inputs, set the count to the maximum count.
 intersect-sum return kmers that occur in all inputs, set the count to the sum of the counts.

 difference return kmers that occur in the first input, but none of the other inputs
 symmetric-difference return kmers that occur in exactly one input

 MODIFIERS:

 output O write kmers generated by the present command to an output meryl database O
 mandatory for count operations.

 EXAMPLES:

 Example: Report 22-mers present in at least one of input1.fasta and input2.fasta.
 Kmers from each input are saved in meryl databases 'input1' and 'input2',
 but the kmers in the union are only reported to the screen.

 meryl print \
 union \
 [count k=22 input1.fasta output input1] \
 [count k=22 input2.fasta output input2]

 Example: Find the highest count of each kmer present in both files, save the kmers to
 database 'maxCount'.

 meryl intersect-max input1 input2 output maxCount

 Example: Find unique kmers common to both files. Brackets are necessary
 on the first 'equal-to' command to prevent the second 'equal-to' from
 being used as an input to the first 'equal-to'.

 meryl intersect [equal-to 1 input1] equal-to 1 input2


```

End the interactive session:

```

[user@cn3316 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





