

document.querySelector('title').textContent = 'Talos: Hyperparameter Optimization for Keras and PyTorch.';
**Talos: Hyperparameter Optimization for Keras and PyTorch.**


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



Talos is a hyperparameter optimization package made for data scientists and data engineers who are tired of mindless parameter hopping and confusing optimization solutions that add complexity instead of reducing it. It works with any Keras, TensorFlow (tf.keras) or PyTorch model, takes minutes to implement, 
involves no new syntax to learn and adds zero new overhead to your workflow. 



### References:


* Autonomio Talos [Computer software]. (2020). Retrieved from http://github.com/autonomio/talos.


Documentation
* [Talos GitHub page](https://github.com/autonomio/talos)


Important Notes
* Module Name: talos (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **TALOS\_HOME**  installation directory
	+ **TALOS\_TEST**  test script directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$  **sinteractive --gres=gpu:k80:1,lscratch=5** 
[user@cn4235 ~]$ **module load talos** 
[+] Loading singularity  3.8.5-1  on cn4235
[+] Loading cuDNN/7.6.5/CUDA-10.1 libraries...
[+] Loading CUDA Toolkit  10.1.105  ...
[+] Loading talos  1.0

```

Copy a test python script, iris.py, to your current foilder:

```

[user@cn4235 ~]$ **cp $TALOS\_TEST/\* .** 

```

Run the test script:

```

[user@cn4235 ~]$**python iris.py** 
...
Successfully opened dynamic library libcublas.so.10
...
 11%|##6                     | 1/9 [00:07<01:02,  7.79s/it]
 22%|#####3                  | 2/9 [00:12<00:42,  6.10s/it]
 33%|########                | 3/9 [00:17<00:33,  5.55s/it]
 44%|##########6             | 4/9 [00:21<00:23,  4.79s/it]
 56%|#############3          | 5/9 [00:24<00:17,  4.37s/it]
 67%|################        | 6/9 [00:28<00:12,  4.13s/it]
 78%|##################6     | 7/9 [00:31<00:07,  3.86s/it]
 89%|#####################3  | 8/9 [00:35<00:03,  3.68s/it]
100%|########################| 9/9 [00:38<00:00,  4.27s/it]

```

End the interactive session:

```

[user@cn4235 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





