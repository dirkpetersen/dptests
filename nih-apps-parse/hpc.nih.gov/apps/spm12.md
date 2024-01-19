

document.querySelector('title').textContent = 'SPM on Biowulf';
SPM on Biowulf


|  |
| --- |
| [SPM12 logo](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)

Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
 |



SPM (Statistical Parametric Mapping) is an interactive application that performs statistical analyses of brain imaging data. 
Statistical Parametric Mapping refers to the construction and assessment of spatially extended statistical processes used to 
test hypotheses about functional imaging data.
The data that can be analyzed by SPM can be a series of images from different cohorts or time-series from the same subject. 
The current release is designed for the analysis of fMRI, PET, SPECT, EEG and MEG. The current version of SPM is SPM12 which has had several updates since 2014.



### References:


* K.J. Friston, A.P. Holmes, K.J. Worsley, J.B. Poline, C. Frith, and R.S.J. Frackowiak. 
 [*Statistical Parametric Maps in Functional Imaging: A General Linear Approach*](https://onlinelibrary.wiley.com/doi/abs/10.1002/hbm.460020402) 
 Human Brain Mapping, 2:189-210, 1995.


Documentation
* [SPM12 main website](https://www.fil.ion.ucl.ac.uk/spm/)
* [SPM12 manual](https://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf)


  

**IMPORTANT: The documentation and examples below refer to the compiled/stand-alone version of SPM.**   

For the uncompiled version of SPM you will need to run Matlab first, then load the SPM toolbox using nih\_matmode load spm. For more details and an example of how to run the uncompiled version of SPM please refer to our [Matlab documentation](/apps/Matlab.html).



Important Notes
* Module Name: spm12 (see [the modules page](/apps/modules.html) for more information)
 * SPM12 is a GUI-based interactive application and requires a graphical connection
to biowulf. Please use the [NoMachine NX client](/docs/nx.html) 
to create the graphical connection as it is the most reliable for SPM12.
 * SPM12 is a compiled MATLAB application so you do not need to have an open matlab session to run it. 
 * It might be necessary to allocate more memory and CPUs than the default shown in the examples below. This depends on the datasets that will be used with SPM12. Please allocate more resources as needed.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn1234 are ready for job

[user@cn1234 ~]$ **module load spm12**
[+] Loading spm12 7870  on cn1234 

[user@cn1234 ~]$ **spm12 &**
Loading SPM12, please wait...

```


You should be able to see the main SPM12 GUI window, and, after choosing a modality, the main SPM12 window for that modality:



![SPM12 Welcome Screen](/images/spm12_GUI.png)


![SPM12 fMRI Screen](/images/spm12_GUI2.png)






