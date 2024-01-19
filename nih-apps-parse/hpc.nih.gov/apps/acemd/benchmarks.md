

document.querySelector('title').textContent = ' Acemd v3 benchmarks';

Acemd v3 benchmarks

**June 2020**

Benchmarks of Acemd v3 test suite on Biowulf GPU nodes. The test suite is provided with the Acemd package, and is available on Biowulf in /usr/local/apps/acemd/acemd3\_examples/.


Hardware:  

* **k80**: Intel(R)\_Xeon(R) CPU E5-2680 v4 @2.40GHz, Tesla K80 
* **p100**: Intel(R)\_Xeon(R) CPU E5-2680 v4 @2.40GHz, Tesla P100-PCIE-16GB 
* **v100**: Intel(R)\_Xeon(R) CPU E5-2680 v4 @2.40GHz, Tesla V100-PCIE-16GB 
* **v100x**: Intel(R) Xeon(R) Gold 6140 CPU @2.30GHz, Tesla V100-SXM2-32GB


Details about the simulation conditions are at [the Acellera website](https://software.acellera.com/docs/latest/acemd3/performance.html). 

The runs were performed with varying numbers of CPUs. In most cases, 2 CPUs per GPU gave the best performance for a given number of GPUs, therefore all results presented below were done with 2 CPUs per GPU.

Based on these benchmarks, there is a significant performance advantage to the v100 and v100x over the p100 and k80 GPUs. However, in most cases there is no advantage to using more than 1 GPU of any type. Before allocating multiple GPUs for your own Acemd jobs, please [run your own benchmarks](/policies/multinode.html#benchmark) to verify that the additional GPUs are providing a performance improvement for your jobs. (and please send any interesting results to staff@hpc.nih.gov).

DHFR (Dihydrofolate reductase)

Number of atoms: 23,558. Box size: 62.23 Å ✕ 62.23 Å ✕ 62.23 Å

![](v3.June2020/dhfr_amber.png)
![](v3.June2020/dhfr_charmm.png)

FactorIX (Factor IX)

Number of atoms: 90,906. Box size: 142.09 Å ✕ 83.34 Å ✕ 78.68 Å

![](v3.June2020/factorIX_amber.png)
![](v3.June2020/factorIX_charmm.png)

STMV (Satellite tobacco mosaic virus)

Number of atoms: 1,067,095. Box size: 221.2 Å ✕ 223.2 Å ✕ 224.5 Å

![](v3.June2020/stmv_amber.png)
![](v3.June2020/stmv_charmm.png)






























