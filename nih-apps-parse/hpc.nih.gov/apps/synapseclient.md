

document.querySelector('title').textContent = 'Synapseclient: an interface to a Synapse ';
**Synapseclient: an interface to a Synapse** 


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



Synapseclient provides an interface to Synapse, a collaborative workspace 
for reproducible data intensive research projects, providing support for:   

(1) integrated presentation of data, code and text,   

(2) fine grained access control, and   

(3) provenance tracking.



Documentation
* [Synapse Home page](https://docs.synapse.org/articles/getting_started.html)
* [Synapse Python Client page](https://python-docs.synapse.org/build/html/index.html)
* [Downloading data](https://docs.synapse.org/articles/downloading_data.html)


Important Notes
* Module Name: synapseclient (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **SYNAPSE\_HOME**  installation directory
	+ **SYNAPSE\_BIN**       executable directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
[user@cn3200 ~]$**module load synapseclient** 
[+] Loading synapseclient  3.0.0
[user@cn3200 ~]$ **synapse -h**
usage: synapse [-h] [--version] [-u SYNAPSEUSER] [-p SYNAPSEPASSWORD]
               [-c CONFIGPATH] [--debug] [-s]
               {get,sync,store,add,mv,cp,associate,delete,query,submit,show,cat,list,set-provenance,get-provenance,set-annotations,get-annotations,create,onweb,login,test-encoding}
               ...

Interfaces with the Synapse repository.

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -u SYNAPSEUSER, --username SYNAPSEUSER
                        Username used to connect to Synapse
  -p SYNAPSEPASSWORD, --password SYNAPSEPASSWORD
                        Password used to connect to Synapse
  -c CONFIGPATH, --configPath CONFIGPATH
                        Path to configuration file used to connect to Synapse
                        [default: /home/user/.synapseConfig]
  --debug
  -s, --skip-checks     suppress checking for version upgrade messages and
                        endpoint redirection

commands:
  The following commands are available:

  {get,sync,store,add,mv,cp,associate,delete,query,submit,show,cat,list,set-provenance,get-provenance,set-annotations,get-annotations,create,onweb,login,test-encoding}
                        For additional help: "synapse  -h"
    get                 downloads a file from Synapse
    sync                Synchronize files described in a manifest to Synapse
    store               uploads and adds a file to Synapse
    add                 uploads and adds a file to Synapse
    mv                  Moves a file/folder in Synapse
    cp                  Copies specific versions of synapse content such as
                        files, folders and projects by recursively copying all
                        sub-content
    associate           Associate local files with the files stored in Synapse
                        so that calls to "synapse get" and "synapse show"
                        don't re-download the files but use the already
                        existing file.
    delete              removes a dataset from Synapse
    query               Performs SQL like queries on Synapse
    submit              submit an entity or a file for evaluation
    show                show metadata for an entity
    cat                 prints a dataset from Synapse
    list                List Synapse entities contained by the given Project
                        or Folder. Note: May not be supported in future
                        versions of the client.
    set-provenance      create provenance records
    get-provenance      show provenance records
    set-annotations     create annotations records
    get-annotations     show annotations records
    create              Creates folders or projects on Synapse
    onweb               opens Synapse website for Entity
    login               login to Synapse and (optionally) cache credentials
    test-encoding       test character encoding to help diagnose problems


```

In order to download data from Synapse, [synapse user account](https://www.synapse.org/#RegisterAccount:0) 
is required.   
  
 A sample command to run synapseclient:

```

[user@cn3200 ~]$ **synapse -u <*your\_synapse\_id*> -p <*your\_synapse\_password*> get syn1899498**
Downloading  [####################]100.00%   7.1kB/7.1kB (26.3MB/s) matrix_100_by_4.tsv Done...
    Downloaded file: matrix_100_by_4.tsv

```

Python-based version of synapseclient is also available:

```

[user@cn3200 ~]$  **python**
Python 3.11.5 | packaged by conda-forge | (main, Aug 27 2023, 03:34:09) [GCC 12.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.

>>> **import synapseclient** 
>>> **syn = synapseclient.login(email = "<*your\_synapse\_id*>", password="<*your\_synapse\_password*>")** 
Welcome, user!
>>> **entity = syn.get('syn1899498')**
Downloading  [####################]100.00%   7.1kB/7.1kB (4.6MB/s) matrix_100_by_4.tsv Done...

```

The previous command stored the data in the local cache. Here is how we can print the metadata from this cache and then to store the data in a local folder, named "syn1899498\_folder" (one folder per data item):

```


>>> **print(entity)**
File: matrix_100_by_4.tsv (syn1899498)
  md5=93c712d2f3096b4d4f078ebddc2068f3
  fileSize=7266
  contentType=text/tab-separated-values
  externalURL=None
  cacheDir=/home/user/.synapseCache/117/30117
  files=['matrix_100_by_4.tsv']
  path=/home/user/.synapseCache/117/30117/matrix_100_by_4.tsv
  synapseStore=True
properties:
  concreteType=org.sagebionetworks.repo.model.FileEntity
  createdBy=377358
  createdOn=2013-05-30T21:58:47.103Z
  dataFileHandleId=30117
  etag=593af68e-7ab2-11e9-98fa-026b0a0ad230
  id=syn1899498
  modifiedBy=377358
  modifiedOn=2013-05-30T21:58:47.103Z
  name=matrix_100_by_4.tsv
  parentId=syn1899495
  versionLabel=1
  versionNumber=1
annotations:
>>> **entity = syn.get('syn1899498', downloadLocation="syn1899498\_folder")**

```

End the interactive session:

```

[user@cn3200 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





