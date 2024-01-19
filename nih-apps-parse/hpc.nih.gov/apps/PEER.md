

document.querySelector('title').textContent = 'PEER: probabilistic estimation of expression residuals ';
**PEER: probabilistic estimation of expression residuals** 


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



PEER stands for "probabilistic estimation of expression residuals". It 
is a collection of Bayesian approaches to infer hidden determinants and 
their effects from gene expression profiles using factor analysis methods.



### References:


* Oliver Stegle, Leopold Parts, Richard Durbin, and John Winn   

*A Bayesian Framework to Account for Complex Non-Genetic Factors 
in Gene Expression Levels Greatly Increases Power in eQTL Studies*   

[PLoS Computational Biology](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000770)  2010, **6**(5), e1000770. doi:10.1371/journal.pcbi.1000770   
* Leopold Parts, Oliver Stegle, John Winn, and Richard Durbin   

*Joint Genetic Analysis of Gene Expression Data with Inferred Cellular Phenotypes*   

[PLoS Genetics](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1001276)  2011, **7**(1), e1001276. doi:10.1371/journal.pgen.1001276


Documentation
* [PEER Tutorial on GitHub](https://github.com/PMBio/peer/wiki/Tutorial)


Important Notes
* Module Name: PEER (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **PEER\_HOME**  PEER installation directory
	+ **PEER\_BIN**       PEER executable directory
	+ **PEER\_DOC**       PEER documentation directory
	+ **PEER\_DATA**    PEER test data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive** 
[user@cn3107 ~]$**module load PEER**
[user@cn3107 ~]$**peertool --help**

USAGE: 

   peertool  [--sigma_off ] [--var\_tol ] [--bound\_tol
 ] [--e\_pb ] [--e\_pa ] [--a\_pb ]
 [--a\_pa ] [-i ] [-n ] [--prior ] [-c
 ] [--var\_file ] -f  [-o ]
 [--has\_header] [--add\_mean] [--no\_a\_out] [--no\_z\_out]
 [--no\_w\_out] [--no\_x\_out] [--no\_res\_out] [--] [--version]
 [-h]
Where: 

 --sigma\_off 
 Variance inactive component

 --var\_tol 
 Variation tolerance

 --bound\_tol 
 Bound tolerance

 --e\_pb 
 Eps node prior parameter b

 --e\_pa 
 Eps node prior parameter a

 --a\_pb 
 Alpha node prior parameter b

 --a\_pa 
 Alpha node prior parameter a

 -i , --n\_iter 
 Number of iterations

 -n , --n\_factors 
 Number of hidden factors

 --prior 
 Factor prior file

 -c , --cov\_file 
 Covariate data file

 --var\_file 
 Expression uncertainty (variance) data file

 -f , --file 
 (required) Expression data file

 -o , --out\_dir 
 Output directory

 --has\_header
 Expression and covariates files have a header

 --add\_mean
 Add a covariate to model mean effect

 --no\_a\_out
 No output of weight precision

 --no\_z\_out
 No output of posterior sparsity prior

 --no\_w\_out
 No output of estimated factor weights

 --no\_x\_out
 No output of estimated factors

 --no\_res\_out
 No output of residual values

 --, --ignore\_rest
 Ignores the rest of the labeled arguments following this flag.

 --version
 Displays version information and exits.

 -h, --help
 Displays usage information and exits.

 Probabilistic estimation of expression residuals (PEER)

[user@cn3107 ~]$**cp $PEER\_DATA/\* .** 
[user@cn3107 ~]$**cd peer/examples/data**

```

Basic usage:

```

[user@cn3107 ~]$**peertool -f expression.csv -n** 
	iteration 0/50
	iteration 1/50
	iteration 2/50
	iteration 3/50
	iteration 4/50
Converged (var(residuals)) after 4 iterations

```

If there are measured experimental variables that may contribute to variability in the data, they can be included in the inference, and specified with the -c flag.

```

[user@cn3107 ~]$**peertool -f expression\_sparse.csv -c covs.csv --out\_dir expression\_sparse\_out**
	iteration 0/50
	iteration 1/50
	iteration 2/50
	iteration 3/50
	iteration 4/50
	iteration 5/50
	iteration 6/50
	iteration 7/50
	iteration 8/50
Converged (var(residuals)) after 8 iterations


```

End the interactive session:

```

[user@cn3107 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





