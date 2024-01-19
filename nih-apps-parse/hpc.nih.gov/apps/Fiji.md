

document.querySelector('title').textContent = 'Fiji: an open-source platform for biological-image analysis on Biowulf';
**Fiji: an open-source platform for
biological-image analysis on Biowulf**


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



Fiji is a distribution of the popular open-source software ImageJ focused on biological-image analysis. Fiji uses modern software engineering practices to combine powerful software libraries with a broad range of scripting languages to enable rapid prototyping of image-processing algorithms. Fiji facilitates the transformation of new algorithms into ImageJ plugins that can be shared with end users through an integrated update system.



### References:


* Johannes Schindelin, Ignacio Arganda-Carreras, Erwin Frise, Verena Kaynig,
Mark Longair, Tobias Pietzsch, Stephan Preibisch, Curtis Rueden, Stephan Saalfeld, Benjamin Schmid, Jean-Yves Tinevez, Daniel James White, Volker Hartenstein,
Kevin Eliceiri, Pavel Tomancak and Albert Cardona,   

*Fiji: an open-source platform for biological-image analysis.*    

Nature Methods, 2012, vol. **9**, N7, p. 676-682.


Documentation
* [Fiji Home page](https://fiji.sc)


Important Notes
* Module Name: Fiji (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Environment variables set
	+ **FIJI\_ROOT**  Fiji installation directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. You should see the Fiji logo appear, and then the Fiji console window.   

In order to use Fiji,
you have to create [graphical X11 connection](https://hpc.nih.gov/docs/connect.html) to Biowulf.   

Both MobaXterm (recommended) and NX will work for Windows users, while XQuartz works well for Mac users.

   

Sample session:



```

[user@biowulf]$ **sinteractive**
[user@@cn3316 ~]$ **module load Fiji** 
[user@@cn3316 ~]$ **fiji &**

```


![Fiji display](/images/Fiji.png)





