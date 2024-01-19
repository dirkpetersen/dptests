

document.querySelector('title').textContent = 'Singularity';


Singularity





|  |
| --- |
| [singularity logo](https://www.sylabs.io/singularity/) |




  


Quick Links
[Important Notes](#notes)
[Creating Singularity containers](#create)
[Binding external directories](#bind)



[Singularity as an Installation Medium](#bind-stationary)


[Interactive Singularity containers](#int)
[Singularity containers in batch](#batch)
[Singularity containers on GPU nodes](#gpu)
[Using Docker containers with Singularity](#docker)
[Troubleshooting containers that hang](#oomkills)



*Extreme Mobility of Compute*

Singularity containers let users run applications in a Linux environment of their choosing. 

Possible uses for Singularity on Biowulf:
* Run an application that was built for a different distribution of Linux than the host OS.
* Reproduce an environment to run a workflow created by someone else.
* Run a series of applications (a 'pipeline') that includes applications built on different platforms.
* Run an application from [Docker Hub](https://hub.docker.com/) on Biowulf without actually installing anything.




**Web sites**
* [Singularity home](https://www.sylabs.io/singularity/)
* [Singularity Documentation](https://sylabs.io/guides/latest/user-guide/)* [Singularity on GitHub](https://github.com/sylabs/singularity)
* [Singularity on Google groups](https://groups.google.com/a/lbl.gov/forum/#!forum/singularity)
* [Docker Hub](https://hub.docker.com/)
* [Singularity Container Services](https://cloud.sylabs.io/home)


**Additional Learning Resources**
* Class taught by NIH HPC staff:
	+ [Materials](https://github.com/NIH-HPC/Singularity-Tutorial/tree/2020-03-10)+ [Day 1 recording](https://youtu.be/3Yg7XI39H4U)+ [Day 2 recording](https://youtu.be/yi82PC--F2U)* [Singularity Basics videos by Sylabs](https://www.youtube.com/playlist?list=PL052H4iYGzysewYEelldGPOgKRJkxd5zp)
* [NIH HPC Singularity Tutorial](https://singularity-tutorial.github.io/)


**Example definition files written by the NIH HPC staff**

These definition files can all be found [on GitHub](https://github.com/NIH-HPC/singularity-examples), and the containers built from them are hosted on [Singularity hub](https://www.singularity-hub.org/collections/267/).
* [**DIGITS** (From the latest NVIDIA Docker Hub image with drivers for NIH HPC GPU nodes)](https://github.com/NIH-HPC/singularity-examples/blob/master/Singularity.digits)* [**Keras/tensorflow**](https://github.com/NIH-HPC/singularity-examples/tree/master/keras)* [**RStudio**](https://github.com/NIH-HPC/singularity-examples/tree/master/rstudio)
* [**Theano** (GPU/CPU support *but requires a GPU to be present either way*, running under Ubuntu 16.04)](https://github.com/NIH-HPC/singularity-examples/blob/master/Singularity.theano)


 Additionally, a large number of staff maintained definition files and associated helper scripts can be found at [this GitHub repo](https://github.com/NIH-HPC/singularity-def-files). These are files that staff members use to install containerized apps on the NIH HPC systems.



Important Notes
[back to top](singularity.html)  












