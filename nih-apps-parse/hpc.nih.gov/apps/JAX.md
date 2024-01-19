

document.querySelector('title').textContent = 'JAX: Autograd and XLA, brought together for high-performance machine learning research. ';
**JAX: Autograd and XLA, brought together for high-performance machine learning research.** 


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



JAX is a python library that brings Autograd and XLA (Accelerated Linear Algebra) together 
for high-performance machine learning research. 
JAX uses XLA to compile and run your NumPy programs on GPUs. 
Compilation happens under the hood by default, with library calls getting just-in-time compiled and executed. But JAX also lets you just-in-time compile your own Python functions into XLA-optimized kernels using a one-function API, jit



Documentation
* [JAX Github page](https://github.com/google/jax)
* [JAX Documentatiom](https://jax.readthedocs.io/en/latest/notebooks/quickstart.html)
* [JAX Tutorial](https://jax.readthedocs.io/en/latest/jax-101/index.html)


Important Notes
* Module Name: JAX (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **JAX\_HOME**  installation directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
  


```

[user@biowulf]$ **sinteractive --gres=gpu:v100:1** 
[user@cn0861 ~]$ **module load jax** 
[+] Loading singularity  3.8.5-1  on cn4207
[+] Loading cuDNN/8.2.1/CUDA-11.3 libraries...
[+] Loading CUDA Toolkit  11.2.2  ...
[+] Loading JAX 0.3.2  ...

```


```

[user@cn0861 ~]$ **python** 

[user@cn0861 ~]$ **>>> import jax.numpy as jnp** 

[user@cn0861 ~]$ **>>> x = jnp.arange(3.0)** 
[user@cn0861 ~]$ **>>> x** 
DeviceArray([0., 1., 2.], dtype=float32)

[user@cn0861 ~]$ **>>> a = jnp.array([0., 1., 3., 4., 5., 4., 3., 2., 1., 0., 0., 0.])** 
[user@cn0861 ~]$ **>>> b = jnp.array([0., 0., 0., 1., 3., 4., 5., 4., 3., 2., 1., 0.])** 
[user@cn0861 ~]$ **>>> jnp.correlate(a, b)** 
DeviceArray([61.], dtype=float32)

[user@cn0861 ~]$ **>>> jnp.correlate(a, b, mode = "same")** 
DeviceArray([26., 43., 61., 75., 81., 75., 61., 43., 26., 13.,  5.,  1.],            dtype=float32)
[user@cn0861 ~]$ **>>> jnp.convolve(a, b, mode = "same")** 
DeviceArray([ 6., 17., 34., 54., 70., 79., 78., 68., 52., 35., 20., 10.],            dtype=float32)


```

etc.


```

[user@cn0861 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226

```





