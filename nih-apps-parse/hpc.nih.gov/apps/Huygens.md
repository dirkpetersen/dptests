

document.querySelector('title').textContent = 'Huygens';
Huygens

![Huygens logo](Huygens_logo.png)
Description
-----------


The Huygens Suite is a complete and easy to use microscopic image restoration (deconvolution) package. Whether applied to images from a wide field (conventional) microscopes, confocal microscopes, 2-photon microscopes, the software in The Huygens Suite will remove blur from your images, improve resolution and reduce noise. It is capable of recovering lost details in very noisy confocal images. Applied to wide field images it will drastically remove blur and improve resolution.


How to Use
----------



Huygens requires an X-windows session. [Click here for more information about X-Windows](/docs/connect.html).


Huygens is licensed to run on a single node in the Biowulf cluster, and only with fixed allocated resources. There
 are three "levels" of allocation: **Small, Medium, and Large**. The necessary level depends on your
 processing needs. Each level allocates a fixed number of GPUs and CPUs. To run Huygens, allocate an
 interactive session with the chosen level, along with an appropriate amount of memory. The memory allocated can be
 any amount up to 240g (GB).


* **Small Allocation** (with 64 GB of memory): 
 
```
[biowulf]$ sinteractive --constraint=huygens --gres=gpu:k80:**1** --cpus-per-task=**14** --mem=64g
```
* **Medium Allocation** (with 128 GB of memory): 
 
```
[biowulf]$ sinteractive --constraint=huygens --gres=gpu:k80:**2** --cpus-per-task=**28** --mem=128g
```
* **Large Allocation** (with 192 GB of memory):
 
```
[biowulf]$ sinteractive --constraint=huygens --gres=gpu:k80:**3** --cpus-per-task=**42** --mem=192g
```


By default, an interactive session is given **8 hours**. If you need more than that amount of time to finish
 your processing, you should include **--time** in your allocation. For example, to allocate 16 hours of
 time to your **Medium** session, type:



```
[biowulf]$ sinteractive --constraint=huygens --gres=gpu:k80:2 --cpus-per-task=28 --mem=128g **--time=16:00:00**
```

Once an interactive session has been started, load the Huygens module:



```
[cn4243]$ module load Huygens
[+] Loading Huygens 19.10 ...
[cn4243]$ huygenspro
```

When finished, make sure to exit the program and the interactive session:



```
[cn4243]$ exit
[biowulf]$
```

There are multiple versions of Huygens available. An easy way of selecting the version is to use modules. To see the modules available, type



```
module avail Huygens
```

To select a module, type



```
module load Huygens/**[ver]**
```

 where **[ver]** is the version of choice. This will set your $PATH variable to allow the Huygens executables.


Type the following commands at the prompt to start the corresponding programs:


* *huygenspro* -- Huygens Professional
* *huscript* -- Huygens Scripting


Disk Space
----------


By default, huygenspro writes to your /home directory. Since /home is limited is disk space, it would be best to create a symlink to your /data directory from the beginning:



```
mkdir /data/**$USER**/SVI
ln -s /data/**$USER**/SVI /home/**$USER**/SVI

```

Documentation
-------------


Click on the 'Help' menu on the far right of the Huygens window after launching the application.


Web-accessible: [SVI](http://www.svi.nl) (Scientific Volume Imaging). Also see the [SVI Support Knowledge Base](http://support.svi.nl).





