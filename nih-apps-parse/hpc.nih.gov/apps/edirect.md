

document.querySelector('title').textContent = 'Entrez Direct E-utilites';
Entrez Direct E-utilites
Description
-----------



Entrez Direct (EDirect) is an advanced method for accessing the NCBI's set of interconnected databases (publication, sequence, structure, gene, variation, expression, etc.) from a UNIX terminal window. Functions take search terms from command-line arguments. Individual operations are combined to build multi-step queries. Record retrieval and formatting normally complete the process.




EDirect also provides an argument-driven function that simplifies the extraction of data from document summaries or other results that are returned in structured XML format. This can eliminate the need for writing custom software to answer ad hoc questions. Queries can move seamlessly between EDirect commands and UNIX utilities or scripts to perform actions that cannot be accomplished entirely within Entrez.



How to Use
----------


Entrez Direct uses [environment modules](modules.html). Type



```
module load edirect
```

at the prompt.


To see the help menu, type



```
efetch -help
```

at the prompt.


**NOTE:** The Entrez Direct utilities will work on Helix, the Biowulf login node, or the Biowulf compute nodes.
Examples
--------


Find the Pubmed entries relating to "Opsin Gene Conversion", then look for all related neighbours of the original papers.


```
esearch -db pubmed -query "opsin gene conversion" | elink -related
```


Find all protein sequences for 'lycopene cyclase' and then fetch those sequences in Fasta format.

```
esearch -db protein -query "lycopene cyclase" | efetch -format fasta
```

And a more complicated example:



```
esearch -db pubmed -query "Beadle AND Tatum AND Neurospora" | \
  elink -related | \
  efilter -query "NOT historical article [FILT]" | \
  efetch -format docsum | \
  xtract -pattern DocumentSummary -present Author -and Title \
    -element Id -first "Author/Name" -element Title | \
  grep -i -e enzyme -e synthesis | \
  sort -t $'\t' -k 2,3f | \
  column -s $'\t' -t | \
  head -n 10 | \
  cut -c 1-80
```

Documentation
-------------


* [Entrez Direct E-utilities Home](http://www.ncbi.nlm.nih.gov/books/NBK179288/)








