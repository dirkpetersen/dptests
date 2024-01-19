

document.querySelector('title').textContent = 'Jupyter on Biowulf';

 .hl { background-color: #ffff99; }
 code {padding: 1px; background-color: #eeeeee; border: 1px solid #bbbbbb;}
 dt {font-weight: bold; margin-top: 5px;}
 dd {padding-left: 2px; border-left: 1px solid #bbbbbb;}

Jupyter on Biowulf


|  |
| --- |
| 
Quick Links
[Pitfalls](#pitfalls)
[Jupyter on biowulf](#general)
[VSCode and Jupyter](https://hpc.nih.gov/apps/vscode.html#jupyter)
 |


Jupyter is an open-source web application that allows you to create and
share documents that contain live code, equations, visualizations and
explanatory text. Uses include: data cleaning and transformation, numerical
simulation, statistical modeling, machine learning and much more. 


### Documentation


* [Jupyter main site](http://jupyter.org/)
* [Video Guide to Jupyter on Biowulf](https://youtu.be/bgLJb1anNPA):
 Includes both a very brief quickstart and more detailed instructions.


Versions
We try to maintain Jupyter kernel configurations for all of our system Python and R modules.
If you require additional kernels or the use of custom environments, you can 
[add your own kernels](https://docs.jupyter.org/en/latest/install/kernels.html)
to your personal configuration.


New versions may add or remove minor additional features or integrations not mentioned on
this page depending on compatibility with newer versions of core packages.


Common pitfalls
This is a list of common pitfalls for jupyter users on the cluster. Some of these
are discussed in other sections of this documentation in more detail.



Running jupyter-notebook on the biowulf login node
It is considered bad form to run jupyter-notebook on the login node. Please
allocate an interactive session to use jupyter-notebook.
Unable to use R magic with python kernels
If the rpy2 module is not loaded before starting the jupyter notebook
 `%load_ext rpy2.ipython`` causes an error. To fix please restart
 jupyter *after* loading the rpy2 module
PermissionError: [Errno 13] Permission denied: '/run/user/xxxx'
This can be fixed by unsetting the `$XDG_RUNTIME_DIR` variable:

```

$ unset XDG_RUNTIME_DIR    # (bash)
$ unsetenv XDG_RUNTIME_DIR # (tcsh)

```

The cause of this error is the `$XDG_RUNTIME_DIR` variable
exported from the login node pointing to a non-existing directory on the 
compute node. This appears to be a problem under only some circumstances.

Expected packages are missing from my python kernel
Jupyter itself is installed as a conda environment and is
the default kernel for new notebooks. This environment, however,
does not include the commonly used scientific stack of packages.
Instead, please use one of the kernels corresponding to the
centrally installed applications (e.g. python/3.10).
Other python modules
No python or R modules need to be loaded. The main python
and R installations are available as pre-configured kernels in
our jupyter setup.
Unable to export PDF
Using the built-in nbconvert to export notebooks to PDF may fail
with the error: "nbconvert failed: xelatex not found on PATH".
To fix this, load the `tex` module.

Starting a Jupyter Instance and connecting from your local browser
**Please see [our Jupyter on Biowulf video](https://www.youtube.com/watch?v=bgLJb1anNPA)
for both a quick start and detailed video guide** to connecting to Jupyter running in a
compute node interactive session. The video shows a Windows client but is equally applicable to macOS or Linux.



In order to connect to a Jupyter Notebook running on a compute node with the browser on
your computer, it is necessary to establish a tunnel from your computer to
biowulf and from biowulf to the compute node. 
[sinteractive](/docs/userguide.html#int) can automatically
create the second leg of this tunnel (from biowulf to the compute node) when started with
the `-T/--tunnel` option. For more details and information on how to set up the
second part of the tunnel see our [tunneling documentation.](/docs/tunneling/)



**Successfully accessing Jupyter on the cluster from your browser is a three step process.
 You must ensure that you do all of the following:**


1. Start your sinteractive session with the `-T`/`--tunnel`
 option. When you start this session up it will give you a specific
 command you can use to access one port on that compute node in the format: `ssh -L
 #####:localhost:##### user@biowulf.nih.gov`
2. Start the Jupyter Lab session with --port $PORT1. This will make sure
 Jupyter is listening on the port that is being forwarded.
3. In a **separate terminal session** from your computer, run the full `ssh
 -L` command from step 1. This is what allows your local browser to actually talk to
 the Jupyter server on the compute node (forwarding that web traffic from
 your local machine, through Biowulf, to the compute node), otherwise,
 you will get a connection refused error.



Note that the python environment hosting the jupyter install is a minimal
python environment. Please use the fully featured kernels named similarly
to their modules (e.g. python/3.10 or R/4.3). The kernels set up their environment
so there is no need to load separate R or python modules before starting jupyter.



Allocate an [interactive session](/docs/userguide.html#int) and
start a jupyter instance as shown below. First, we launch [tmux](https://github.com/tmux/tmux/wiki) (or [screen](https://www.gnu.org/software/screen/)) on the login node so
that we don't lose our session if our connection to the login node drops.

```

[user@biowulf]$ **module load tmux** *# You can use screen instead; you don't need to module load it*
[user@biowulf]$ **tmux**
[user@biowulf]$ **sinteractive --gres=lscratch:5 --mem=10g --tunnel**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job
                                                                           
Created 1 generic SSH tunnel(s) from this compute node to                  
biowulf for your use at port numbers defined                               
in the $PORTn ($PORT1, ...) environment variables.                         
                                                                           
                                                                           
Please create a SSH tunnel from your workstation to these ports on biowulf.
On Linux/MacOS, open a terminal and run:                                   
                                                                           
    ssh  -L 33327:localhost:33327 biowulf.nih.gov                          
                                                                           
For Windows instructions, see https://hpc.nih.gov/docs/tunneling           

[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load jupyter**
[user@cn3144]$ # if you want to use R magick with python kernels you also need the rpy2 module
[user@cn3144]$ **module load rpy2**
[user@cn3144]$ **cp ${JUPYTER\_TEST\_DATA:-none}/\* .**
[user@cn3144]$ **ls -lh**
total 196K
-rw-r--r-- 1 user group 7.8K Oct  1 14:31 Pokemon.csv
-rw-r--r-- 1 user group 186K Oct  1 14:31 Seaborn_test.ipynb

[user@cn3144]$ **jupyter kernelspec list**
Available kernels:
  bash                     /usr/local/apps/jupyter/conda/envs/5.3.0/share/jupyter/kernels/bash
  ir42                     /usr/local/apps/jupyter/conda/envs/5.3.0/share/jupyter/kernels/ir42
  ir43                     /usr/local/apps/jupyter/conda/envs/5.3.0/share/jupyter/kernels/ir43
  jupyter_matlab_kernel    /usr/local/apps/jupyter/conda/envs/5.3.0/share/jupyter/kernels/jupyter_matlab_kernel
  py3.10                   /usr/local/apps/jupyter/conda/envs/5.3.0/share/jupyter/kernels/py3.10
  py3.8                    /usr/local/apps/jupyter/conda/envs/5.3.0/share/jupyter/kernels/py3.8
  py3.9                    /usr/local/apps/jupyter/conda/envs/5.3.0/share/jupyter/kernels/py3.9
  python3                  /usr/local/apps/jupyter/conda/envs/5.3.0/share/jupyter/kernels/python3
  sos                      /usr/local/apps/jupyter/conda/envs/5.3.0/share/jupyter/kernels/sos

```

In the example show here I will use the convenient environment variable `$PORT1` which 
expands to the port reserved by sinteractive - 33327 in this example. You can also use the port number
directly. Brief explanation of the options:



`notebook|lab|console`
Start a jupyter notebook, jupyter lab, or jupyter console
`--ip localhost`
listen on localhost only
`--port $PORT1`
listen on `$PORT1`, the variable set when using `--tunnel` for sinteractive
`--no-browser`
We're not running the browser on the compute node!

Option 1: Running Jupyter Lab


Jupyter Lab is an improved interface with additional capabilities and features.



```

[user@cn3144]$ **jupyter lab --ip localhost --port $PORT1 --no-browser** 
[I 14:06:30.925 ServerApp] ipyparallel | extension was successfully linked.
[I 14:06:30.931 ServerApp] jupyter_server_mathjax | extension was successfully linked.
[I 14:06:30.931 ServerApp] jupyter_server_proxy | extension was successfully linked.
[I 14:06:30.940 ServerApp] jupyterlab | extension was successfully linked.
[I 14:06:30.940 ServerApp] jupyterlab_git | extension was successfully linked
...
[C 2022-11-21 14:06:31.112 ServerApp]

    To access the server, open this file in a browser:
        file:///spin1/home/linux/user/.local/share/jupyter/runtime/jpserver-62322-open.html
    Or copy and paste one of these URLs:
        http://localhost:40792/lab?token=xxxxxxxxxx
     or http://127.0.0.1:40792/lab?token=xxxxxxxxxx

```

 Keep this open for as long as you're using your notebook. **At this point you must make sure
 you have a tunnel from your PC to this specific Jupyter session.**


Option 2: Running Jupyter Notebook


Jupyter Notebook is the simple, classic web interface to Jupyter.



```

[user@cn3144]$ **jupyter notebook --ip localhost --port $PORT1 --no-browser** 
[I 12:48:25.645 NotebookApp] [nb_conda_kernels] enabled, 20 kernels found
[I 12:48:26.053 NotebookApp] [nb_anacondacloud] enabled
[I 12:48:26.077 NotebookApp] [nb_conda] enabled
[I 12:48:26.322 NotebookApp] ✓ nbpresent HTML export ENABLED
[W 12:48:26.323 NotebookApp] ✗ nbpresent PDF export DISABLED: No module named nbbrowserpdf.exporters.pdf
[I 12:48:26.330 NotebookApp] Serving notebooks from local directory: /spin1/users/user
[I 12:48:26.330 NotebookApp] 0 active kernels 
[I 12:48:26.330 NotebookApp] The Jupyter Notebook is running at: http://localhost:33327/?token=xxxxxxxxxx
[I 12:48:26.331 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[C 12:48:26.333 NotebookApp] 
    
    Copy/paste this URL into your browser when you connect for the first time,
    to login with a token:
        http://localhost:33327/?token=xxxxxxxxxx


```


For documentation on how to connect a browser on your local computer to your jupyter instance on
a biowulf compute node via a tunnel see our general
[tunneling documentation](/docs/tunneling/).



### Selecting a kernel


**In the notebook interface** the kernels highlighted in red below correspond
to the equally named modules on the command line interface. Note that no modules other than
jupyter has to be loaded - the kernels do all the required setup:



![Select notebook kernels](/images/jupyter_notebook_kernels.png)

Similarly, **in the Jupyter Lab interface**



![Select lab kernels](/images/jupyter_lab_kernels.png)

Kernels may change as new python or R installations are added or old ones are retired. Look
for kernels with the same name as command line modules.






