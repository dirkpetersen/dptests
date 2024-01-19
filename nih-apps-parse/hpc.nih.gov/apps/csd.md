

document.querySelector('title').textContent = 'Cambridge Structural Database/WebCSD';
Cambridge Structural Database/WebCSD
![CCDC logo](/images/csd.jpg)

The Cambridge Structural Database is the world repository of small molecule crystal structures. The database holds bibliographic, 2D chemical and 3D structural results from crystallographic analyses of organics, organometallics and metal complexes. Both X-ray and neutron diffraction studies are included for compounds containing up to ca. 500 atoms (including hydrogens). NOTE: The database is NOW available on the biowulf cluster.


Along with the database are tools for search, retreival, analysis, and display of the CSD contents. These include ConQuest (for text, numeric, 2D substructure and 3D geometric searching) and VISTA (statistical analysis and display of geometric and other data). 





|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Web CSD](#web)
[Interactive job](#int) 
[Python API - interactive](#pint)
[Python API - batch job](#pbatch) 
 |


Installers for Windows, Mac and Linux OSes, are also available for all NIH employees interested in using the software. If you would like to use CSD from your local desktop, please contact staff@hpc.nih.gov, for a link to the download.


Documentation

* [CSD and Conquest documentation](http://www.ccdc.cam.ac.uk/support/documentation/#csds) at the CSD website.
* [What is WebCSD](http://www.ccdc.cam.ac.uk/products/csd_system/webcsd/)* [Interactve WebCSD Demo](http://webcsd.ccdc.cam.ac.uk/teaching_database_demo.php)




WebCSD

NIH users can access the CSD directly via the [WebCSD portal](https://www.ccdc.cam.ac.uk/structures/). New structures are added weekly for up-to-date searches. Users can perform searches - Similarity, Substructure, Text and Numerical - as well as browse the database.

The portal is restricted to NIH IP addresses only. If you are at home, you must connect to the NIH VPN before using.


Interactive Quest 

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


on Biowulf
CSD requires an X-windows, or preferably NX, session. [Click here for more information about X-Windows and NX](/docs/connect.html).


To use, type *module load CSD*, then *quest&*, at the prompt.



Sample session: (replace 'username' with your own Helix/Biowulf username).

You cannot run Conquest on the Biowulf login node. You will first need to allocate an interactive session with 'sinteractive'.


```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load CSD**

[user@cn3144 ~]$  **quest &**
Starting ConQuest. Use questv5 for QUEST.
Running ConQuest with glibc version 2.5
2.5

```


![](/images/conquest1.jpg)

![](/images/Conquest2.png)










Python API on Biowulf

Loading the CSD module will also set up the paths for the CSD python API installation. Sample interactive session:

```

biowulf% **sinteractive --mem=20g** 
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load CSD**
[+] Loading Cambridge Structural Database 2020  on cn0931

# Note: the following is required to run without an Xwindows display
# but will prevent diagram generation and other graphics-oriented
# pieces of the CSD from running
[user@cn3144 ~]$ **export CCDC\_PYTHON\_API\_NO\_QAPPLICATION=True**


[user@cn3144 ~]$ **python**
Python 3.7.6 | packaged by conda-forge | (default, Jan  7 2020, 22:33:48)
[GCC 7.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> **from ccdc import io**
QStandardPaths: XDG_RUNTIME_DIR not set, defaulting to '/tmp/runtime-user'
>>> **csd\_reader = io.EntryReader('CSD')**
>>> c**ryst\_abebuf = csd\_reader.crystal('ABEBUF')**
>>> **mol\_abebuf = csd\_reader.molecule('ABEBUF')**
>>> **round(mol\_abebuf.molecular\_weight, 3)**
317.341
>>> **mol\_abebuf.is\_organic**
True
>>> **print(mol\_abebuf.heaviest\_component.smiles)**
O=C1Nc2ccccc2C(=O)Nc2ccccc12
>>> **quit()**

[user@cn3144 ~]$ **exit**
biowulf% 

```



Python API in batch job

Set up a batch script along the following lines:

```

#!/bin/bash

module load CSD

# Note: the following is required to run without an Xwindows display
# but will prevent diagram generation and other graphics-oriented
# pieces of the CSD from running
export CCDC_PYTHON_API_NO_QAPPLICATION=True

python << EOF

from ccdc import io
from ccdc.search import TextNumericSearch
csd_reader = io.EntryReader('CSD')
text_numeric_search = TextNumericSearch()
text_numeric_search.add_compound_name('aspirin')
identifiers = [h.identifier for h in text_numeric_search.search()]
for identifier in sorted(set(identifiers)):	
      e = csd_reader.entry(identifier)
      if e.melting_point:
	    print('%-8s http://dx.doi.org/%-25s %s' % (e.identifier,e.publication.doi,e.melting_point))

EOF

```


Submit with:

```

sbatch myscript

```

The script should produce an output file containing the following:

```

[+] Loading Cambridge Structural Database 2020  on cn####
ACMEBZ   http://dx.doi.org/10.1107/S0567740881006729 385K
ACSALA13 http://dx.doi.org/10.1021/ja056455b         135.5deg.C
BEHWOA   http://dx.doi.org/10.1107/S0567740882005731 401K
CUASPR01 http://dx.doi.org/10.1107/S1600536803026126 above 573 K
HUPPOX   http://dx.doi.org/10.1039/b208574g          91-96 deg.C
HUPPOX01 http://dx.doi.org/None                      91-96 deg.C
NUWTOP01 http://dx.doi.org/10.1039/c2ce06313a        147 deg.C
PIKYUG   http://dx.doi.org/10.1021/acs.cgd.8b00718   502 K

```



















