

document.querySelector('title').textContent = 'Apptainer (the Linux Foundation variant of Singularity)';


Apptainer (the Linux Foundation variant of Singularity)





|  |
| --- |
| [Apptainer logo](https://apptainer.org/) |




  

Quick Links

[Documentation](#doc)
[Important Notes](#notes)
[Interactive Apptainer containers](#int)
[Binding external directories](#bind)
[Containers on GPUs](#gpu)
[Containers in batch](#batch)
[Swarms of containers](#swarm)
[Troubleshooting containers](#oomkills)
[Creating containers](#create)
[Faking app installations](#apps)
[Extended example](#example)





Apptainer is a tool allowing you to build and run [Linux containers](https://opensource.com/resources/what-are-linux-containers). Linux containers can be thought of as small, lightweight virtual machines encapsulating an entire operating system. Containers let users run applications in a Linux environment of their choosing.

Possible uses for Apptainer on Biowulf:
* Run a containerized application on Biowulf without actually installing anything. Prebuilt containers can be obtained from a variety of locations such as:
	+ [Docker Hub](https://hub.docker.com/)+ [Quay.io](https://quay.io/search)+ [NVIDIA NGC](https://catalog.ngc.nvidia.com/containers)+ [Sylabs Cloud Library](https://cloud.sylabs.io/library)+ [BioContainers](https://biocontainers.pro/registry)* Run an application that was built for a different distribution of Linux than the host OS.
* Reproduce an environment to run a workflow created by someone else.
* Run a series of applications (a 'pipeline') that includes applications built on different platforms.



Please note, Apptainer gives you the ability to install and run applications in your own Linux environment with your own customized software stack. With this ability comes the added responsibility of managing your own Linux environment. While the NIH HPC staff can provide guidance on how to create and use Apptainer containers, we do not have the resources to manage containers for individual users. If you decide to use Apptainer, it is your responsibility to manage your own containers.


  

  


Documentation
[back to top](apptainer.html)  










