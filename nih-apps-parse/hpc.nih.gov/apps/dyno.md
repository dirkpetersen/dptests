

document.querySelector('title').textContent = 'dyno: inferring single-cell trajectories';
**dyno: inferring single-cell trajectories**


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



dyno is a meta package that installs several other packages
from the [dynverse](https://github.com/dynverse).   

It comprises a set of R packages 
to construct and interpret single-cell trajectories. 



### References:


* Wouter Saelens, Robrecht Cannoodt, Helena Todorov & Yvan Saeys  

*A comparison of single-cell trajectory inference methods*    

[Nature Biotechnology](https://www.nature.com/articles/s41587-019-0071-9) **37**, 
547-554 (2019).


Documentation
* [dyno Github page](https://github.com/dynverse/dyno)


Important Notes
* Module Name: dyno (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Environment variables set
	+ **DYNOHOME**  dyno installation directory
	+ **R\_LIBS**   dyno R library directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf ~]$ **sinteractive --mem=12g -c4 --gres=lscratch:10 --tunnel**
...
Created 1 generic SSH tunnel(s) from this compute node to
biowulf for your use at port numbers defined
in the $PORTn ($PORT1, ...) environment variables.

Please create a SSH tunnel from your workstation to these ports on biowulf.
On Linux/MacOS, open a terminal and run:

    ssh  -L 33521:localhost:33521 user@biowulf.nih.gov

For Windows instructions, see https://hpc.nih.gov/docs/tunneling
[user@cn0868 ~]$

```

Store the port number you get (in this case: 33521) and the id of the compute node you logged in to, in this case cn0868  
 
On your local system (PC or Mac), open a (second) terminal/linux shell and type:

```
 
**ssh -t -L 33521:localhost:33521 biowulf "ssh -L 33521:localhost:33521 cn0868"** 

```

(make sure to use the actual port number you stored instead of 33521 and the actual node id number instead of cn0868).   

Now return to the original (first) terminal/linux shell.   


```
 
[user@cn0868 ~]$ **module load dyno**
[+] Loading dyno  20220906  on cn0868
[+] Loading singularity  3.8.0  on cn0868

[user@cn0868 ~]$ **R**
>  **library(dyno)** 
Loading required package: dynfeature
Loading required package: dynguidelines
Loading required package: dynmethods
Loading required package: dynplot
Loading required package: dynwrap
>  **data("fibroblast\_reprogramming\_treutlein")**  
>  **dataset <- wrap\_expression(counts = fibroblast\_reprogramming\_treutlein$counts, expression = fibroblast\_reprogramming\_treutlein$expression)** 
>  **guidelines <- guidelines\_shiny(dataset, port=33521, launch.browser=F)** 
Loading required package: shiny

Listening on http://0.0.0.0:33521
Warning: All elements of `...` must be named.
Did you want `renderers = c(column_id, renderer, label, title, style, default, name, trajectory_type,
    category_old, scaling_type)`?
Loading required namespace: akima


```
 
On your local system (PC or Mac), open a web browser or new tab on it, and navigate to:   


```

 **http://localhost:33521** 

```

(make sure to use the port number you stored instead of 33521)   
   

![](dyno/dynguidelines.png)
  
   

Select method(s) to be used, then close the browser window/tab, or click on "Close & use"

```

>  **methods\_selected <- guidelines$methods\_selected** 
>  **methods\_selected** 
[1] "slingshot" "paga_tree" "scorpius"  "angle"
>  **model <- infer\_trajectory(dataset, methods\_selected[1])** 
Running singularity exec 'docker://dynverse/ti_slingshot:v1.0.3' echo hi
...
Running /usr/local/current/singularity/3.8.0/bin/singularity exec \
  --containall -B \
  '/tmp/RtmppFshCi/file941a3d4b0be4/:/copy_mount,/tmp/RtmppFshCi/file941a8f3b018/tmp:/tmp2' \
  'docker://dynverse/ti_slingshot:v1.0.3' cp /code/definition.yml \
  /copy_mount/
...

```

(the computation may take a few minutes)

```

>  **library(magrittr)** 
>  **model <- model %>% add\_dimred(dyndimred::dimred\_mds, expression\_source = dataset$expression)**
>  **plot\_dimred(model, expression\_source = dataset$expression, grouping = fibroblast\_reprogramming\_treutlein$grouping )** 
Coloring by grouping
Loading required namespace: RColorBrewer

```

![](dyno/plot_dimread.png)

img {
 display: block;
 margin-left: auto;
 margin-right: auto;
}

  


```

>  **q()** 
[user@cn0868 ~]$ **exit**
exit
salloc.exe: Relinquishing job allocation 49998864

[user@biowulf ~]$ 

```





