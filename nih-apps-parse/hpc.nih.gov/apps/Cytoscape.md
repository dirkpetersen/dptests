

document.querySelector('title').textContent = 'Cytoscape on Biowulf';
Cytoscape on Biowulf
[![Cytoscape Logo](/images/cy3logoOrange.svg)](http://www.cytoscape.org)

 Cytoscape is an open source software platform for visualizing complex networks and integrating
 these with any type of attribute data. With Cytoscape, you can load molecular and genetic interaction
 data sets in many standards formats, project and integrate global datasets and functional annotations,
 establish powerful visual mappings across these data, perform advanced analysis and modeling and
 visualize and analyze human-curated pathway datasets.





|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int)
[Cytoscape Configuration](#config)
[Cytoscape Apps](#apps)
 |


Documentation
* [Cytoscape Home](http://www.cytoscape.org)* [Cytoscape User Manual](http://manual.cytoscape.org/en/stable/index.html)


Important Notes
This application requires a [graphical connection using NX](/docs/nx.html)


* Module Name: Cytoscape (see [the modules page](/apps/modules.html) for more information)
* Environment variables set: CYTOSCAPE\_HOME
* Example files in $CYTOSCAPE\_HOME/sampleData/
* Refrain from manually loading java modules as these may cause problems with Cytoscape
* We recommend allocated at least 8GB memory for Cytoscape jobs and sessions
* To use the cyREST API with client applications like the RCy3 library, make sure to set the following before running Cytoscape:
 export no\_proxy=localhost



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) from Biowulf using a graphical connection (we recommend NX) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=6**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

```

To use Cytoscape, you must load the module. Type: 




```

**module load Cytoscape**
**Cytoscape&**

```

![Cytoscape GUI](/images/cy3_3_0_desktop.png)


Cytoscape Configuration Directory
Cytoscape will create a configuration directory inside the user's home directory. This is usually 
 /home/$USER/CytoscapeConfiguration. We have observed two potential issues around this:


* Several batch jobs that automate task(s) in Cytoscape, usually running via an API, can write many files to
 the configuration directory. In some cases this can exceed the 16GB quota on the user's home directory causing
 many problems to the user and running jobs. A possible solution is to reset JAVA\_OPTIONS and change the
 configuration directory location. For example, run the following before starting Cytoscape:



```
export JAVA_OPTIONS="-Duser.home=/data/$USER/CytoscapeConfiguration"
```
* If user's switch between different versions of the Cytoscape module, they may encounter some incompatibilities
 since the same configuration directory is used regardless of which Cytoscape version is started. This could be
 especially true for Cytoscape plugins (apps) that are incompatible with certain versions of Cytoscape. Use the
 same solution as listed above or rename the previous configuration directory.


Cytoscape Apps
Cytoscape functionality can be extended using plugins, called Cytoscape Apps. There is an 
 "[app store](https://apps.cytoscape.org/)" from where users can download, install, and manage these
 apps. However, to use the app store from a compute node on Biowulf (i.e. from an interactive session) Cytoscape
 must be configured to use the node's [http proxy](https://hpc.nih.gov/docs/transfer.html#compute).
 Here are the steps on how to do this:



1. Allocate an sinteractive session as usual and run the following **before starting Cytoscape**:



```
echo $http_proxy
```


An example of the output you should see after running the above command is: http://dtn05-e0:3128.
 Here the proxy server is *dtn05-e0* and port is *3128*. You will likely see different values. Note them down.

- Run Cytoscape and navigate to **Edit** --> **Preferences**> --> **Proxy Settings**. For "Type" select **http**. Enter the
 proxy server and port you obtained in the previous steps. Leave the User Name and Password fields blank.
 Click "OK".

- Then open Cytoscape's App Manager by navigating to **Apps** --> **App Manager**. You should see the
 app manager contact the default App Store, https://apps.cytoscape.org and refresh the cache of
 available apps for installation.



 If you want to install and manage Cytoscape apps in a different sinteractive session, you will may need to redo step 1.
 This is because the proxy server may be different between different compute nodes.





