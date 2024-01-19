

document.querySelector('title').textContent = 'LJA: Assembling Long and Accurate Reads Using Multiplex de Bruijn Graphs ';
**LJA: Assembling Long and Accurate Reads Using Multiplex de Bruijn Graphs** 


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



The La Jolla Assembler (LJA) is a tool for assemblies of long and accurate reads. It reduces the error rate in
these reads by three orders of magnitude (making them nearly error-free) and constructs the de Bruijn
graph for large genomes and large k-mer sizes. Since the de Bruijn graph constructed for a fixed k-mer
size is typically either too tangled or too fragmented, LJA uses a new concept of a multiplex de Bruijn
graph with varying k-mer sizes. 



### References:


* Anton Bankevich, Andrey Bzikadze, Mikhail Kolmogorov, Dmitry Antipov, Pavel A. Pevzner   

*Multiplex de Bruijn graphs enable genome assembly from long, high-fidelity reads*   

[Nat Biotechnol (2022). https://doi.org/10.1038/s41587-022-01220-6](https://www.nature.com/articles/s41587-022-01220-6)


Documentation
* [LJA on GitHub](https://github.com/AntonBankevich/LJA)


Important Notes
* Module Name: LJA (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **LJA\_HOME**  installation directory
	+ **LJA\_BIN**     executable directory
	+ **LJA\_SRC**     source code directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive** 
[user@cn3107 ~]$**module load LJA**

```

Available executables:

```

[user@biowulf]$  **ls $LJA\_BIN** 
dot_bulge_stats  jumboDBG  lja  run_polishing  run_tests  sdbg_stats

```

Basic usage:

```

[user@cn3107 ~]$**lja -h**
LJA: genome assembler for PacBio HiFi reads based on de Bruijn graph.
Usage: lja [options] -o <output-dir> --reads <reads_file> [--reads <reads_file2> ...]

Basic options:
  -o <file_name> (or --output-dir <file_name>)  Name of output folder. Resulting graph will be stored there.
  --reads <file_name>                           Name of file that contains reads in fasta or fastq format. This option can be used any number of times in the same command line. In this case reads from all specified files will be used as an input.
  -h (or --help)                                Print this help message.

Advanced options:
  -t <int> (or --threads <int>)                 Number of threads. The default value is 16.
  -k <int>                                      Value of k used for initial error correction.
  -K <int>                                      Value of k used for final error correction and initialization of multiDBG.
  --diploid                                     Use this option for diploid genomes. By default LJA assumes that the genome is haploid or inbred.

```

Running tests:

```

[user@cn3107 ~]$ **run\_tests** 
Running main() from /usr/local/apps/LJA/LJA/src/tests/googletest-master/googletest/src/gtest_main.cc
[==========] Running 30 tests from 26 test suites.
[----------] Global test environment set-up.
[----------] 1 test from DB1
[ RUN      ] DB1.Basic
[       OK ] DB1.Basic (0 ms)
[----------] 1 test from DB1 (0 ms total)

[----------] 1 test from DBSingleEdge1
[ RUN      ] DBSingleEdge1.Basic
[       OK ] DBSingleEdge1.Basic (0 ms)
[----------] 1 test from DBSingleEdge1 (0 ms total)

[----------] 1 test from DBSingleEdge2
[ RUN      ] DBSingleEdge2.Basic
[       OK ] DBSingleEdge2.Basic (0 ms)
[----------] 1 test from DBSingleEdge2 (0 ms total)

[----------] 1 test from DBSingleEdge3
[ RUN      ] DBSingleEdge3.Basic
[       OK ] DBSingleEdge3.Basic (0 ms)
[----------] 1 test from DBSingleEdge3 (0 ms total)

[----------] 1 test from DBStVertex
[ RUN      ] DBStVertex.Basic
[       OK ] DBStVertex.Basic (0 ms)
[----------] 1 test from DBStVertex (0 ms total)

[----------] 1 test from DBEvVertex
[ RUN      ] DBEvVertex.Basic
[       OK ] DBEvVertex.Basic (0 ms)
[----------] 1 test from DBEvVertex (0 ms total)

[----------] 2 tests from DB1inVertex
[ RUN      ] DB1inVertex.Basic
[       OK ] DB1inVertex.Basic (1 ms)
[ RUN      ] DB1inVertex.WithShortEdge
[       OK ] DB1inVertex.WithShortEdge (0 ms)
[----------] 2 tests from DB1inVertex (1 ms total)

[----------] 2 tests from DB1outVertex
[ RUN      ] DB1outVertex.Basic
[       OK ] DB1outVertex.Basic (0 ms)
[ RUN      ] DB1outVertex.WithShortEdge
[       OK ] DB1outVertex.WithShortEdge (0 ms)
[----------] 2 tests from DB1outVertex (0 ms total)

[----------] 1 test from DBComplexVertex
[ RUN      ] DBComplexVertex.Basic
[       OK ] DBComplexVertex.Basic (0 ms)
[----------] 1 test from DBComplexVertex (0 ms total)

[----------] 1 test from DBComplexVertexLoop1
[ RUN      ] DBComplexVertexLoop1.Basic
[       OK ] DBComplexVertexLoop1.Basic (0 ms)
[----------] 1 test from DBComplexVertexLoop1 (0 ms total)

[----------] 1 test from DBComplexVertexLoop2
[ RUN      ] DBComplexVertexLoop2.Basic
[       OK ] DBComplexVertexLoop2.Basic (0 ms)
[----------] 1 test from DBComplexVertexLoop2 (0 ms total)

[----------] 1 test from DBComplexVertexLoop3
[ RUN      ] DBComplexVertexLoop3.Basic
[       OK ] DBComplexVertexLoop3.Basic (0 ms)
[----------] 1 test from DBComplexVertexLoop3 (0 ms total)

[----------] 1 test from DBComplexVertexLoop4
[ RUN      ] DBComplexVertexLoop4.Basic
[       OK ] DBComplexVertexLoop4.Basic (0 ms)
[----------] 1 test from DBComplexVertexLoop4 (0 ms total)

[----------] 1 test from DBComplexVertexLoop5
[ RUN      ] DBComplexVertexLoop5.Basic
[       OK ] DBComplexVertexLoop5.Basic (1 ms)
[----------] 1 test from DBComplexVertexLoop5 (1 ms total)

[----------] 1 test from DBComplexVertexConn4
[ RUN      ] DBComplexVertexConn4.Basic
[       OK ] DBComplexVertexConn4.Basic (0 ms)
[----------] 1 test from DBComplexVertexConn4 (0 ms total)

[----------] 1 test from DBComplexVertexConn3
[ RUN      ] DBComplexVertexConn3.Basic
[       OK ] DBComplexVertexConn3.Basic (0 ms)
[----------] 1 test from DBComplexVertexConn3 (0 ms total)

[----------] 1 test from DBComplexVertexConn3_2
[ RUN      ] DBComplexVertexConn3_2.Basic
[       OK ] DBComplexVertexConn3_2.Basic (0 ms)
[----------] 1 test from DBComplexVertexConn3_2 (0 ms total)

[----------] 1 test from DBComplexVertexLoop6
[ RUN      ] DBComplexVertexLoop6.Basic
[       OK ] DBComplexVertexLoop6.Basic (0 ms)
[----------] 1 test from DBComplexVertexLoop6 (0 ms total)

[----------] 1 test from DBComplexVertexLoop7
[ RUN      ] DBComplexVertexLoop7.Basic
[       OK ] DBComplexVertexLoop7.Basic (0 ms)
[----------] 1 test from DBComplexVertexLoop7 (0 ms total)

[----------] 1 test from DBComplexVertexLoop8
[ RUN      ] DBComplexVertexLoop8.Basic
[       OK ] DBComplexVertexLoop8.Basic (0 ms)
[----------] 1 test from DBComplexVertexLoop8 (0 ms total)

[----------] 1 test from DBDoubleLoopRC
[ RUN      ] DBDoubleLoopRC.Basic
[       OK ] DBDoubleLoopRC.Basic (1 ms)
[----------] 1 test from DBDoubleLoopRC (1 ms total)

[----------] 1 test from DBIsolate
[ RUN      ] DBIsolate.Basic
[       OK ] DBIsolate.Basic (0 ms)
[----------] 1 test from DBIsolate (0 ms total)

[----------] 1 test from DBEmptyGraph
[ RUN      ] DBEmptyGraph.Basic
[       OK ] DBEmptyGraph.Basic (0 ms)
[----------] 1 test from DBEmptyGraph (0 ms total)

[----------] 2 tests from RRPathsTest
[ RUN      ] RRPathsTest.Basic
[       OK ] RRPathsTest.Basic (0 ms)
[ RUN      ] RRPathsTest.MergeIterDereference
[       OK ] RRPathsTest.MergeIterDereference (0 ms)
[----------] 2 tests from RRPathsTest (0 ms total)

[----------] 1 test from EdgeSegment
[ RUN      ] EdgeSegment.Basic
[       OK ] EdgeSegment.Basic (0 ms)
[----------] 1 test from EdgeSegment (0 ms total)

[----------] 2 tests from MDBGSeq
[ RUN      ] MDBGSeq.SingleSegment
[       OK ] MDBGSeq.SingleSegment (0 ms)
[ RUN      ] MDBGSeq.SeveralSegments
[       OK ] MDBGSeq.SeveralSegments (0 ms)
[----------] 2 tests from MDBGSeq (0 ms total)

[----------] Global test environment tear-down
[==========] 30 tests from 26 test suites ran. (3 ms total)
[  PASSED  ] 30 tests.

```

End the interactive session:

```

[user@cn3107 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





