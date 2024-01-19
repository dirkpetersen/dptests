/data/

document.querySelector('title').textContent = ' Amber22 benchmarks';

Amber22 benchmarks

**April 2023**



|  |
| --- |
| 
Quick Links
[CPU benchmarks](#cpu)
[GPU benchmarks](#gpu)
 |


The benchmark runs below used the Amber 20 Benchmark suite, downloadable from [here](http://ambermd.org/GPUPerformance.php). 

Older benchmarks: [[Amber 20](amber20.html)] [[Amber 18](amber18B.html)] [Amber 16](amber16.html)] [[Amber 14](amber14.html)] 

Based on these benchmarks, there is a significant performance advantage to running Amber on the GPU nodes rather than on a CPU-only node.


CPU-only Benchmarks

Amber 22 with all patches as of July 2022, built with gcc 9.2.0, OpenMPI 4.1.3 (Biowulf module 'amber/22.gcc'), running under Rocky8   

Hardware: 64 x 2.8 GHz (AMD Epyc 7543)

**Explicit Solvent**  

[![](22/explicit.cpus.png)](22/explicit.cpus.png)  

**Implicit Solvent**  

[![](22/myoglobin.cpus.png)](22/myoglobin.cpus.png)
[![](22/nucleosome.cpus.png)](22/nucleosome.cpus.png)



The benchmark above was run on a AMD Epyc 7543 @2.8GHz with 64 cores (128 hyperthreaded cores). As with other 
MD applications, the performance drops when more than 64 MPI processes 
(i.e. 1 process per physical core) are run. Thus, if running on CPUs, it is important to use the '--ntasks-per-core=1' flag 
when submitting the job, to ensure that only 1 MPI process is run on each physical core. If it is possible to run on a GPU node, 
you will get significantly better performance as shown in the benchmarks below.


GPU benchmarks

Amber 22 with all patches as of July 2022, built with gcc 9.2.0, CUDA 11.3, and OpenMPI 4.1.3 (Biowulf module 'amber/22-gpu'), running under Rocky8  

Hardware:
* k80: Intel(R)\_Xeon(R)\_CPU\_E5-2680\_v4\_@\_2.40GHz, Tesla K80 -- benchmark runs with 1-4 GPUs on single node
* p100: Intel(R)\_Xeon(R)\_CPU\_E5-2680\_v4\_@\_2.40GHz, Tesla\_P100-PCIE-16GB -- benchmark runs with 1-4 GPUs on single node
* v100: Intel(R)\_Xeon(R)\_CPU\_E5-2680\_v4\_@\_2.40GHz, Tesla\_V100-PCIE-16GB -- benchmark runs with 1-4 GPUs on single node
* v100x: Intel(R)\_Xeon(R)\_Gold\_6140@\_2.30GHz, Tesla\_V100-SXM2-32GB -- benchmark runs with 1-4 GPUs on single node
* a100: AMD EPYC 7543P@\_2.40GHz, Tesla\_A100-SXM4-80GB -- benchmark runs with 1-4 GPUs on single node



**Explicit Solvent (PME)**
![](22/Cellulose_NPT_4fs.png)
![](22/Cellulose_NVE_4fs.png)
![](22/FactorIX_NPT_4fs.png)
![](22/FactorIX_NVE_4fs.png)
![](22/JAC_NPT_4fs.png)
![](22/JAC_NVE_4fs.png)
![](22/STMV_NPT_4fs.png)
![](22/STMV_NVE_4fs.png)

**Implicit Solvent (GB)**
![](22/myoglobin.png)
![](22/nucleosome.png)
![](22/TRPCage.png)


























