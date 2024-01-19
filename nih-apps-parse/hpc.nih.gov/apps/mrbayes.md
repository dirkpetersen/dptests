

document.querySelector('title').textContent = 'MrBayes on Biowulf';
MrBayes on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 

 |



 MrBayes is a program for Bayesian inference and model choice across a wide range of phylogenetic and evolutionary models. MrBayes uses Markov chain Monte Carlo (MCMC) methods to estimate the posterior distribution of model parameters.



### References:


* Ronquist, F., M. Teslenko, P. van der Mark, D.L. Ayres, A. Darling, S. Höhna, B. Larget, L. Liu, M.A. Suchard, and J.P. Huelsenbeck. 2012. MRBAYES 3.2: Efficient Bayesian phylogenetic inference and model selection across a large model space. Syst. Biol. 61:539-542. [doi:10.1093/sysbio/sys029](https://doi.org/10.1093/sysbio/sys029)
* Ronquist, F., and J.P. Huelsenbeck. 2003. MRBAYES 3: Bayesian phylogenetic inference under mixed models. Bioinformatics 19:1572-1574. [doi:10.1093/bioinformatics/btg180](https://doi.org/10.1093/bioinformatics/btg180)
* Huelsenbeck, J.P., and F. Ronquist. 2001. MRBAYES: Bayesian inference of phylogeny. Bioinformatics 17:754-755. [doi:10.1093/bioinformatics/17.8.754](https://doi.org/10.1093/bioinformatics/17.8.754)


Documentation
* [MrBayes Main Site](http://nbisweden.github.io/MrBayes)
* [Manual](http://nbisweden.github.io/MrBayes/manual.html)
* [MrBayes on GitHub](https://github.com/NBISweden/MrBayes)


Important Notes
* Module Name: mrbayes (see [the modules page](/apps/modules.html) for more information)
* Can run single-threaded, [multi-threaded with MPI](/docs/userguide.html#par), or on a single [GPU](/docs/userguide.html#gpu). From the manual:
 
 The GPU code can be a lot faster than the CPU code, particularly for amino
 acid and codon models. However, the length of the sequences also influences the
 speed-up. In general, the longer the sequences are, the better GPU performance
 you can expect. If the sequences are short, the overhead involved in shuffling
 data to and from the GPU may well overshadow any performance gain you get in
 the computation step. Try the various calculator options out in short runs before
 you decide on the best option for longer runs.
 



 If you run multiple simultaneous runs or use Metropolis coupling, which is
 standard practice with MrBayes, then you can speed up your analyses
 considerably by using the MPI version of MrBayes.
* MrBayes supports **checkpointing**, so use sbatch's --no-requeue flag to allow editing your job to resume from a checkpoint in case of failure.
 From the manual:
 
 Since version 3.2, MrBayes prints all parameter values of all chains (cold and
 heated) to a checkpoint file every Checkfreq generations, by default every 2,000
 generations. The checkpoint file has the suffix .ckp. If you run an analysis and it
 is stopped prematurely, you can restart it from the last checkpoint by using mcmc append=yes.
 MrBayes will start the new analysis from the checkpoint; it will
 even read in all the old trees and include them in the convergence diagnostics. At
 the end of the new run, you will have parameter and tree files that are
 indistinguishable from those you would have obtained from an uninterrupted
 analysis.
 

* Environment variables set 
	+ MRBAYES\_HOME* Example files in $MRBAYES\_HOME/share/examples/mrbayes/



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.

 The following sample sessions run the simple tutorial described in the manual using single-threaded, multi-threaded, and GPU modes.
 They differ in how the allocation is created and, in the case of MPI operation, how the application is first started.





 single CPU
 




```

[user@biowulf user]$ **sinteractive**
salloc.exe: Pending job allocation 24780751
salloc.exe: job 24780751 queued and waiting for resources
salloc.exe: job 24780751 has been allocated resources
salloc.exe: Granted job allocation 24780751
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3160 are ready for job
srun: error: x11: no local DISPLAY defined, skipping
[user@cn3160 user]$ **module load mrbayes**
[+] Loading gcc  6.4.0  ... 
[+] Loading openmpi 3.1.2  for GCC 6.4.0 
[+] Loading mrbayes, version 3.2.7... 
[user@cn3160 user]$ **mkdir mrbayes && cd mrbayes**
[user@cn3160 mrbayes]$ **cp $MRBAYES\_HOME/share/examples/mrbayes/primates.nex .**
[user@cn3160 mrbayes]$ **mb**


                            MrBayes 3.2.7 x86_64

                      (Bayesian Analysis of Phylogeny)

                             (Parallel version)
                         (1 processors available)

              Distributed under the GNU General Public License


               Type "help" or "help <command>" for information
                     on the commands that are available.

                   Type "about" for authorship and general
                       information about the program.


MrBayes > **execute primates.nex**

   Executing file "primates.nex"
   DOS line termination
   Longest line length = 920
   Parsing file
   Expecting NEXUS formatted file
   Reading data block
      Allocated taxon set
      Allocated matrix
      Defining new matrix with 12 taxa and 898 characters
      Data is Dna
      Data matrix is not interleaved
      Gaps coded as -
      Taxon  1 -> Tarsius_syrichta
      Taxon  2 -> Lemur_catta
      Taxon  3 -> Homo_sapiens
      Taxon  4 -> Pan
      Taxon  5 -> Gorilla
      Taxon  6 -> Pongo
      Taxon  7 -> Hylobates
      Taxon  8 -> Macaca_fuscata
      Taxon  9 -> M_mulatta
      Taxon 10 -> M_fascicularis
      Taxon 11 -> M_sylvanus
      Taxon 12 -> Saimiri_sciureus
      Successfully read matrix
      Setting default partition (does not divide up characters)
      Setting model defaults
      Seed (for generating default start values) = 1555622841
      Setting output file names to "primates.nex.run<i>.<p|t>"
   Exiting data block
   Reached end of file

MrBayes > **lset nst=6 rates=invgamma**

   Setting Nst to 6
   Setting Rates to Invgamma
   Successfully set likelihood model parameters

MrBayes > **mcmc ngen=20000 samplefreq=100**
printfreq=100 diagnfreq=1000
   Setting number of generations to 20000
   Setting sample frequency to 100
   Running Markov chain
   MCMC stamp = 8692445572
   Seed = 30626700
   Swapseed = 1555622841
   Model settings:

      Data not partitioned --
         Datatype  = DNA
         Nucmodel  = 4by4
         Nst       = 6
                     Substitution rates, expressed as proportions
                     of the rate sum, have a Dirichlet prior
                    (1.00,1.00,1.00,1.00,1.00,1.00)
         Covarion  = No
         # States  = 4
                     State frequencies have a Dirichlet prior
                     (1.00,1.00,1.00,1.00)
         Rates     = Invgamma
                     The distribution is approximated using 4 categories.
                     Shape parameter is exponentially
                     distributed with parameter (1.00).
                     Proportion of invariable sites is uniformly dist-
                     ributed on the interval (0.00,1.00).

   Active parameters: 

      Parameters
      ---------------------
      Revmat              1
      Statefreq           2
      Shape               3
      Pinvar              4
      Ratemultiplier      5
      Topology            6
      Brlens              7
      ---------------------

      1 --  Parameter  = Revmat
            Type       = Rates of reversible rate matrix
            Prior      = Dirichlet(1.00,1.00,1.00,1.00,1.00,1.00)

      2 --  Parameter  = Pi
            Type       = Stationary state frequencies
            Prior      = Dirichlet

      3 --  Parameter  = Alpha
            Type       = Shape of scaled gamma distribution of site rates
            Prior      = Exponential(1.00)

      4 --  Parameter  = Pinvar
            Type       = Proportion of invariable sites
            Prior      = Uniform(0.00,1.00)

      5 --  Parameter  = Ratemultiplier
            Type       = Partition-specific rate multiplier
            Prior      = Fixed(1.0)

      6 --  Parameter  = Tau
            Type       = Topology
            Prior      = All topologies equally probable a priori
            Subparam.  = V

      7 --  Parameter  = V
            Type       = Branch lengths
            Prior      = Unconstrained:GammaDir(1.0,0.1000,1.0,1.0)


   Number of chains per MPI processor = 8

   The MCMC sampler will use the following moves:
      With prob.  Chain will use move
          0.93 %   Dirichlet(Revmat)
          0.93 %   Slider(Revmat)
          0.93 %   Dirichlet(Pi)
          0.93 %   Slider(Pi)
          1.85 %   Multiplier(Alpha)
          1.85 %   Slider(Pinvar)
          9.26 %   ExtSPR(Tau,V)
          9.26 %   ExtTBR(Tau,V)
          9.26 %   NNI(Tau,V)
          9.26 %   ParsSPR(Tau,V)
         37.04 %   Multiplier(V)
         12.96 %   Nodeslider(V)
          5.56 %   TLMultiplier(V)

	    Division 1 has 413 unique site patterns
	    Initializing conditional likelihoods

	    Running benchmarks to automatically select fastest BEAGLE resource... 
	    Using BEAGLE v3.1.2 resource 0 for division 1:
	    Rsrc Name : CPU
	    Impl Name : CPU-4State-Single
	    Flags: PROCESSOR_CPU PRECISION_SINGLE COMPUTATION_SYNCH EIGEN_REAL
            SCALING_MANUAL SCALERS_RAW VECTOR_NONE THREADING_CPP      MODEL STATES: 4
	    Initializing invariable-site conditional likelihoods

	    Initial log likelihoods and log prior probs for run 1:
	    Chain 1 -- -8803.174992 -- 42.620562
	    Chain 2 -- -8887.607821 -- 42.620562
	    Chain 3 -- -8612.869815 -- 42.620562
	    Chain 4 -- -8522.921441 -- 42.620562

	    Initial log likelihoods and log prior probs for run 2:
	    Chain 1 -- -8679.106083 -- 42.620562
	    Chain 2 -- -9046.716149 -- 42.620562
	    Chain 3 -- -8486.154983 -- 42.620562
	    Chain 4 -- -8760.764335 -- 42.620562


	    Using a relative burnin of 25.0 % for diagnostics

	    Chain results (20000 generations requested):

	    0 -- [-8803.175] (-8887.608) (-8612.870) (-8522.921) * [-8679.106] (-9046.716) (-8486.155) (-8760.764) (...0 remote chains...) 
	    1000 -- (-5992.439) (-6049.599) (-6013.840) [-5968.341] * [-5911.315] (-6029.790) (-6070.369) (-6060.764) (...0 remote chains...) -- 0:00:19
	    2000 -- [-5864.131] (-5894.405) (-5927.265) (-5891.980) * [-5785.441] (-5916.646) (-6005.412) (-5875.297) (...0 remote chains...) -- 0:00:36
	    3000 -- (-5819.796) (-5839.428) (-5889.297) [-5787.887] * [-5741.932] (-5828.153) (-5866.805) (-5775.116) (...0 remote chains...) -- 0:01:25
	    4000 -- [-5741.953] (-5797.157) (-5773.872) (-5751.430) * (-5743.271) (-5806.691) (-5820.159) [-5744.328] (...0 remote chains...) -- 0:01:04
	    5000 -- [-5729.026] (-5781.991) (-5748.914) (-5735.459) * (-5731.240) (-5750.312) (-5766.323) [-5734.272] (...0 remote chains...) -- 0:00:51

	    Average standard deviation of split frequencies: 0.000000

	    6000 -- [-5726.299] (-5758.040) (-5738.510) (-5724.524) * (-5722.869) [-5738.720] (-5757.223) (-5723.481) (...0 remote chains...) -- 0:00:42
	    7000 -- (-5724.209) (-5726.639) [-5726.048] (-5723.930) * (-5722.860) (-5723.057) (-5741.001) [-5725.379] (...0 remote chains...) -- 0:01:12
	    8000 -- (-5722.836) (-5738.059) [-5726.183] (-5735.521) * (-5721.778) (-5727.580) (-5735.357) [-5726.899] (...0 remote chains...) -- 0:01:00
	    9000 -- (-5722.199) [-5724.381] (-5729.216) (-5727.787) * (-5724.959) (-5728.901) (-5737.999) [-5729.190] (...0 remote chains...) -- 0:00:50
	    10000 -- (-5722.838) (-5724.077) [-5723.598] (-5723.395) * (-5728.865) [-5729.452] (-5722.015) (-5733.628) (...0 remote chains...) -- 0:00:42

	    Average standard deviation of split frequencies: 0.001034

	    11000 -- [-5722.571] (-5742.559) (-5728.971) (-5726.905) * (-5730.962) (-5727.398) (-5724.377) [-5729.311] (...0 remote chains...) -- 0:00:34
	    12000 -- (-5728.377) [-5722.619] (-5724.693) (-5723.195) * (-5723.553) [-5719.051] (-5734.158) (-5729.540) (...0 remote chains...) -- 0:00:28
	    13000 -- (-5728.154) (-5730.298) (-5722.790) [-5726.406] * (-5729.549) [-5721.158] (-5739.788) (-5741.450) (...0 remote chains...) -- 0:00:23
	    14000 -- (-5725.711) [-5728.830] (-5725.949) (-5723.317) * (-5731.914) (-5724.208) (-5730.943) [-5729.157] (...0 remote chains...) -- 0:00:19
	    15000 -- [-5724.016] (-5743.916) (-5729.716) (-5726.448) * (-5725.390) (-5729.683) (-5728.894) [-5720.638] (...0 remote chains...) -- 0:00:15

	    Average standard deviation of split frequencies: 0.001378

	    16000 -- (-5726.544) (-5732.291) (-5726.437) [-5729.115] * (-5729.030) (-5731.877) (-5722.591) [-5719.019] (...0 remote chains...) -- 0:00:12
	    17000 -- (-5723.574) [-5724.367] (-5721.819) (-5726.479) * (-5728.811) (-5728.894) [-5723.818] (-5733.027) (...0 remote chains...) -- 0:00:08
	    18000 -- (-5735.050) (-5723.310) [-5725.305] (-5726.638) * (-5723.073) (-5728.121) (-5720.994) [-5731.460] (...0 remote chains...) -- 0:00:05
	    19000 -- [-5731.532] (-5725.784) (-5723.220) (-5716.554) * (-5740.380) (-5718.239) (-5730.634) [-5727.171] (...0 remote chains...) -- 0:00:02
	    20000 -- (-5725.996) (-5726.782) (-5720.546) [-5721.486] * [-5721.475] (-5724.620) (-5726.671) (-5717.234) (...0 remote chains...) -- 0:00:00

	    Average standard deviation of split frequencies: 0.001041

Continue with analysis? (yes/no): **no**

            Analysis completed in 52 seconds
	    Analysis used 25.60 seconds of CPU time on processor 0
	    Likelihood of best state for "cold" chain of run 1 was -5715.13
	    Likelihood of best state for "cold" chain of run 2 was -5715.39

	    Acceptance rates for the moves in the "cold" chain of run 1:
	    With prob.   (last 100)   chain accepted proposals by move
            36.9 %     ( 40 %)     Dirichlet(Revmat)
            61.2 %     ( 61 %)     Slider(Revmat)
            17.3 %     ( 19 %)     Dirichlet(Pi)
            24.3 %     ( 27 %)     Slider(Pi)
            37.6 %     ( 34 %)     Multiplier(Alpha)
            65.9 %     ( 69 %)     Slider(Pinvar)
            0.1 %     (  0 %)     ExtSPR(Tau,V)
            0.2 %     (  0 %)     ExtTBR(Tau,V)
            0.2 %     (  0 %)     NNI(Tau,V)
            0.4 %     (  0 %)     ParsSPR(Tau,V)
            37.2 %     ( 25 %)     Multiplier(V)
            23.2 %     ( 22 %)     Nodeslider(V)
            14.7 %     ( 11 %)     TLMultiplier(V)

	    Acceptance rates for the moves in the "cold" chain of run 2:
	    With prob.   (last 100)   chain accepted proposals by move
            33.3 %     ( 36 %)     Dirichlet(Revmat)
            65.3 %     ( 63 %)     Slider(Revmat)
            16.7 %     ( 10 %)     Dirichlet(Pi)
            22.9 %     ( 19 %)     Slider(Pi)
            40.6 %     ( 39 %)     Multiplier(Alpha)
            62.7 %     ( 65 %)     Slider(Pinvar)
            0.3 %     (  0 %)     ExtSPR(Tau,V)
            0.3 %     (  0 %)     ExtTBR(Tau,V)
            0.2 %     (  0 %)     NNI(Tau,V)
            0.5 %     (  1 %)     ParsSPR(Tau,V)
            36.8 %     ( 29 %)     Multiplier(V)
            23.9 %     ( 25 %)     Nodeslider(V)
            15.2 %     ( 16 %)     TLMultiplier(V)

	    Chain swap information for run 1:

            1     2     3     4 
            --------------------------
	    1 |        0.67  0.42  0.29 
	    2 |  3381        0.65  0.45 
	    3 |  3280  3347        0.66 
	    4 |  3323  3332  3337       

	    Chain swap information for run 2:

            1     2     3     4 
            --------------------------
	    1 |        0.66  0.41  0.23 
	    2 |  3374        0.61  0.40 
	    3 |  3268  3287        0.62 
	    4 |  3364  3340  3367       

	    Upper diagonal: Proportion of successful state exchanges between chains
	    Lower diagonal: Number of attempted state exchanges between chains

	    Chain information:

	    ID -- Heat 
	    -----------
	    1 -- 1.00  (cold chain)
	    2 -- 0.91 
	    3 -- 0.83 
	    4 -- 0.77 

	    Heat = 1 / (1 + T * (ID - 1))
	    (where T = 0.10 is the temperature and ID is the chain number)


MrBayes > **quit**

   Quitting program

[user@cn3160 mrbayes]$ **exit**
salloc.exe: Relinquishing job allocation 24780751
[user@biowulf user]$ 

```





 multiple CPUs with MPI
 




```

[user@biowulf ~]$ **sinteractive --ntasks 8 --ntasks-per-core 1**
salloc.exe: Pending job allocation 24837664
salloc.exe: job 24837664 queued and waiting for resources
salloc.exe: job 24837664 has been allocated resources
salloc.exe: Granted job allocation 24837664
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3112 are ready for job
srun: error: x11: no local DISPLAY defined, skipping
[user@cn3112 ~]$ **module load mrbayes**
[+] Loading gcc  6.4.0  ... 
[+] Loading openmpi 3.1.2  for GCC 6.4.0 
[+] Loading mrbayes, version 3.2.7... 
[user@cn3112 ~]$ **mkdir mrbayes && cd mrbayes**
[user@cn3112 mrbayes]$ **cp $MRBAYES\_HOME/share/examples/mrbayes/primates.nex .**
[user@cn3112 mrbayes]$ **mpirun -np $SLURM\_NTASKS mb**


                            MrBayes 3.2.7 x86_64

                      (Bayesian Analysis of Phylogeny)

                             (Parallel version)
                         (8 processors available)

              Distributed under the GNU General Public License


               Type "help" or "help <command>" for information
                     on the commands that are available.

                   Type "about" for authorship and general
                       information about the program.


MrBayes >  **execute primates.nex**

   Executing file "primates.nex"
   DOS line termination
   Longest line length = 920
   Parsing file
   Expecting NEXUS formatted file
   Reading data block
      Allocated taxon set
      Allocated matrix
      Defining new matrix with 12 taxa and 898 characters
      Data is Dna
      Data matrix is not interleaved
      Gaps coded as -
      Taxon  1 -> Tarsius_syrichta
      Taxon  2 -> Lemur_catta
      Taxon  3 -> Homo_sapiens
      Taxon  4 -> Pan
      Taxon  5 -> Gorilla
      Taxon  6 -> Pongo
      Taxon  7 -> Hylobates
      Taxon  8 -> Macaca_fuscata
      Taxon  9 -> M_mulatta
      Taxon 10 -> M_fascicularis
      Taxon 11 -> M_sylvanus
      Taxon 12 -> Saimiri_sciureus
      Successfully read matrix
      Setting default partition (does not divide up characters)
      Setting model defaults
      Seed (for generating default start values) = 1555705864
      Setting output file names to "primates.nex.run<i>.<p|t>"
   Exiting data block
   Reached end of file

MrBayes > **lset nst=6 rates=invgamma**

	    Setting Nst to 6
	    Setting Rates to Invgamma
	    Successfully set likelihood model parameters

MrBayes > **mcmc ngen=20000 samplefreq=100**

	    Setting number of generations to 20000
	    Setting sample frequency to 100
	    Running Markov chain
	    MCMC stamp = 5367313892
	    Seed = 3019911
	    Swapseed = 1555705864
	    Model settings:

	    Data not partitioned --
            Datatype  = DNA
            Nucmodel  = 4by4
            Nst       = 6
            Substitution rates, expressed as proportions
            of the rate sum, have a Dirichlet prior
            (1.00,1.00,1.00,1.00,1.00,1.00)
            Covarion  = No
            # States  = 4
            State frequencies have a Dirichlet prior
            (1.00,1.00,1.00,1.00)
            Rates     = Invgamma
            The distribution is approximated using 4 categories.
            Shape parameter is exponentially
            distributed with parameter (1.00).
            Proportion of invariable sites is uniformly dist-
            ributed on the interval (0.00,1.00).

	    Active parameters: 

	    Parameters
	    ---------------------
	    Revmat              1
	    Statefreq           2
	    Shape               3
	    Pinvar              4
	    Ratemultiplier      5
	    Topology            6
	    Brlens              7
	    ---------------------

	    1 --  Parameter  = Revmat
            Type       = Rates of reversible rate matrix
            Prior      = Dirichlet(1.00,1.00,1.00,1.00,1.00,1.00)

	    2 --  Parameter  = Pi
            Type       = Stationary state frequencies
            Prior      = Dirichlet

	    3 --  Parameter  = Alpha
            Type       = Shape of scaled gamma distribution of site rates
            Prior      = Exponential(1.00)

	    4 --  Parameter  = Pinvar
            Type       = Proportion of invariable sites
            Prior      = Uniform(0.00,1.00)

	    5 --  Parameter  = Ratemultiplier
            Type       = Partition-specific rate multiplier
            Prior      = Fixed(1.0)

	    6 --  Parameter  = Tau
            Type       = Topology
            Prior      = All topologies equally probable a priori
            Subparam.  = V

	    7 --  Parameter  = V
            Type       = Branch lengths
            Prior      = Unconstrained:GammaDir(1.0,0.1000,1.0,1.0)


	    Number of chains per MPI processor = 1

	    The MCMC sampler will use the following moves:
	    With prob.  Chain will use move
            0.93 %   Dirichlet(Revmat)
            0.93 %   Slider(Revmat)
            0.93 %   Dirichlet(Pi)
            0.93 %   Slider(Pi)
            1.85 %   Multiplier(Alpha)
            1.85 %   Slider(Pinvar)
            9.26 %   ExtSPR(Tau,V)
            9.26 %   ExtTBR(Tau,V)
            9.26 %   NNI(Tau,V)
            9.26 %   ParsSPR(Tau,V)
            37.04 %   Multiplier(V)
            12.96 %   Nodeslider(V)
            5.56 %   TLMultiplier(V)

	    Division 1 has 413 unique site patterns
	    Initializing conditional likelihoods

	    Running benchmarks to automatically select fastest BEAGLE resource... 
	    Using BEAGLE v3.1.2 resource 0 for division 1:
	    Rsrc Name : CPU
	    Impl Name : CPU-4State-Single
	    Flags: PROCESSOR_CPU PRECISION_SINGLE COMPUTATION_SYNCH EIGEN_REAL
            SCALING_MANUAL SCALERS_RAW VECTOR_NONE THREADING_CPP      MODEL STATES: 4
	    Initializing invariable-site conditional likelihoods

	    Initial log likelihoods and log prior probs for run 1:
	    Chain 1 -- -8886.343679 -- 42.620562

	    There are 7 more chains on other processor(s)


	    Using a relative burnin of 25.0 % for diagnostics

	    Chain results (20000 generations requested):

	    0 -- [-8886.344] [...7 remote chains...] 
	    1000 -- [-5901.771] [...7 remote chains...] -- 0:00:00
	    2000 -- [-5836.424] [...7 remote chains...] -- 0:00:00
	    3000 -- (-5812.036) [...7 remote chains...] -- 0:00:05
	    4000 -- (-5752.144) [...7 remote chains...] -- 0:00:04
	    5000 -- (-5732.309) [...7 remote chains...] -- 0:00:03

	    Average standard deviation of split frequencies: 0.002015

	    6000 -- [-5722.850] [...7 remote chains...] -- 0:00:04
	    7000 -- (-5728.707) [...7 remote chains...] -- 0:00:03
	    8000 -- (-5723.999) [...7 remote chains...] -- 0:00:03
	    9000 -- (-5729.417) [...7 remote chains...] -- 0:00:03
	    10000 -- (-5728.500) [...7 remote chains...] -- 0:00:02

	    Average standard deviation of split frequencies: 0.003101

	    11000 -- (-5724.098) [...7 remote chains...] -- 0:00:02
	    12000 -- (-5724.540) [...7 remote chains...] -- 0:00:02
	    13000 -- (-5735.209) [...7 remote chains...] -- 0:00:02
	    14000 -- (-5731.311) [...7 remote chains...] -- 0:00:01
	    15000 -- [-5719.576] [...7 remote chains...] -- 0:00:01

	    Average standard deviation of split frequencies: 0.001378

	    16000 -- (-5728.118) [...7 remote chains...] -- 0:00:01
	    17000 -- (-5738.108) [...7 remote chains...] -- 0:00:00
	    18000 -- (-5736.811) [...7 remote chains...] -- 0:00:00
	    19000 -- (-5733.949) [...7 remote chains...] -- 0:00:00
	    20000 -- (-5728.526) [...7 remote chains...] -- 0:00:00

	    Average standard deviation of split frequencies: 0.001561

   Continue with analysis? (yes/no): **no**

	    Analysis completed in 6 seconds
	    Analysis used 6.41 seconds of CPU time on processor 0
	    Likelihood of best state for "cold" chain of run 1 was -5716.44
	    Likelihood of best state for "cold" chain of run 2 was -5716.44

	    Acceptance rates for the moves in the "cold" chain of run 1:
	    With prob.   (last 100)   chain accepted proposals by move
            32.3 %     ( 28 %)     Dirichlet(Revmat)
            61.6 %     ( 61 %)     Slider(Revmat)
            18.1 %     ( 16 %)     Dirichlet(Pi)
            22.1 %     ( 22 %)     Slider(Pi)
            41.1 %     ( 47 %)     Multiplier(Alpha)
            62.9 %     ( 66 %)     Slider(Pinvar)
            0.1 %     (  0 %)     ExtSPR(Tau,V)
            0.2 %     (  0 %)     ExtTBR(Tau,V)
            0.4 %     (  0 %)     NNI(Tau,V)
            0.6 %     (  0 %)     ParsSPR(Tau,V)
            36.8 %     ( 26 %)     Multiplier(V)
            23.8 %     ( 24 %)     Nodeslider(V)
            14.5 %     ( 17 %)     TLMultiplier(V)

	    Acceptance rates for the moves in the "cold" chain of run 2:
	    With prob.   (last 100)   chain accepted proposals by move
            31.5 %     ( 36 %)     Dirichlet(Revmat)
            60.8 %     ( 60 %)     Slider(Revmat)
            19.7 %     ( 24 %)     Dirichlet(Pi)
            21.4 %     ( 21 %)     Slider(Pi)
            33.2 %     ( 29 %)     Multiplier(Alpha)
            71.1 %     ( 71 %)     Slider(Pinvar)
            0.3 %     (  0 %)     ExtSPR(Tau,V)
            0.1 %     (  0 %)     ExtTBR(Tau,V)
            0.6 %     (  3 %)     NNI(Tau,V)
            0.6 %     (  0 %)     ParsSPR(Tau,V)
            36.2 %     ( 40 %)     Multiplier(V)
            24.3 %     ( 21 %)     Nodeslider(V)
            15.3 %     ( 19 %)     TLMultiplier(V)

	    Chain swap information for run 1:

            1     2     3     4 
            --------------------------
	    1 |        0.68  0.46  0.34 
	    2 |  3307        0.69  0.47 
	    3 |  3354  3333        0.65 
	    4 |  3258  3360  3388       

	    Chain swap information for run 2:

            1     2     3     4 
            --------------------------
	    1 |        0.60  0.44  0.32 
	    2 |  3297        0.69  0.48 
	    3 |  3278  3255        0.71 
	    4 |  3334  3353  3483       

	    Upper diagonal: Proportion of successful state exchanges between chains
	    Lower diagonal: Number of attempted state exchanges between chains

	    Chain information:

	    ID -- Heat 
	    -----------
	    1 -- 1.00  (cold chain)
	    2 -- 0.91 
	    3 -- 0.83 
	    4 -- 0.77 

	    Heat = 1 / (1 + T * (ID - 1))
	    (where T = 0.10 is the temperature and ID is the chain number)


MrBayes > **quit**

   Deleting previously defined characters
   Deleting previously defined taxa
   Quitting program

[user@cn3112 mrbayes]$ **exit**
salloc.exe: Relinquishing job allocation 24837664
[user@biowulf ~]$ 
	
```





 single GPU
 




```

[user@biowulf ~]$ **sinteractive --gres gpu:p100:1**
salloc.exe: Pending job allocation 24838742
salloc.exe: job 24838742 queued and waiting for resources
salloc.exe: job 24838742 has been allocated resources
salloc.exe: Granted job allocation 24838742
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn2351 are ready for job
srun: error: x11: no local DISPLAY defined, skipping
[user@cn2351 ~]$ **module load mrbayes**
[+] Loading gcc  6.4.0  ... 
[+] Loading openmpi 3.1.2  for GCC 6.4.0 
[+] Loading mrbayes, version 3.2.7... 
[user@cn2351 ~]$ **mkdir mrbayes && cd mrbayes**
[user@cn2351 mrbayes]$ **cp $MRBAYES\_HOME/share/examples/mrbayes/primates.nex .**
[user@cn2351 mrbayes]$ **mb**


                            MrBayes 3.2.7 x86_64

                      (Bayesian Analysis of Phylogeny)

                             (Parallel version)
                         (1 processors available)

              Distributed under the GNU General Public License


               Type "help" or "help <command>" for information
                     on the commands that are available.

                   Type "about" for authorship and general
                       information about the program.


MrBayes > **execute primates.nex**

   Executing file "primates.nex"
   DOS line termination
   Longest line length = 920
   Parsing file
   Expecting NEXUS formatted file
   Reading data block
      Allocated taxon set
      Allocated matrix
      Defining new matrix with 12 taxa and 898 characters
      Data is Dna
      Data matrix is not interleaved
      Gaps coded as -
      Taxon  1 -> Tarsius_syrichta
      Taxon  2 -> Lemur_catta
      Taxon  3 -> Homo_sapiens
      Taxon  4 -> Pan
      Taxon  5 -> Gorilla
      Taxon  6 -> Pongo
      Taxon  7 -> Hylobates
      Taxon  8 -> Macaca_fuscata
      Taxon  9 -> M_mulatta
      Taxon 10 -> M_fascicularis
      Taxon 11 -> M_sylvanus
      Taxon 12 -> Saimiri_sciureus
      Successfully read matrix
      Setting default partition (does not divide up characters)
      Setting model defaults
      Seed (for generating default start values) = 1555707139
      Setting output file names to "primates.nex.run<i>.<p|t>"
   Exiting data block
   Reached end of file

MrBayes > **lset nst=6 rates=invgamma**

   Setting Nst to 6
   Setting Rates to Invgamma
   Successfully set likelihood model parameters

MrBayes > **mcmc ngen=20000 samplefreq=100**

   Setting number of generations to 20000
   Setting sample frequency to 100
   Running Markov chain
   MCMC stamp = 5005632409
   Seed = 1921241950
   Swapseed = 1555707139
   Model settings:

      Data not partitioned --
         Datatype  = DNA
         Nucmodel  = 4by4
         Nst       = 6
                     Substitution rates, expressed as proportions
                     of the rate sum, have a Dirichlet prior
                     (1.00,1.00,1.00,1.00,1.00,1.00)
         Covarion  = No
         # States  = 4
                     State frequencies have a Dirichlet prior
                     (1.00,1.00,1.00,1.00)
         Rates     = Invgamma
                     The distribution is approximated using 4 categories.
                     Shape parameter is exponentially
                     distributed with parameter (1.00).
                     Proportion of invariable sites is uniformly dist-
                     ributed on the interval (0.00,1.00).

   Active parameters: 

      Parameters
      ---------------------
      Revmat              1
      Statefreq           2
      Shape               3
      Pinvar              4
      Ratemultiplier      5
      Topology            6
      Brlens              7
      ---------------------

      1 --  Parameter  = Revmat
            Type       = Rates of reversible rate matrix
            Prior      = Dirichlet(1.00,1.00,1.00,1.00,1.00,1.00)

      2 --  Parameter  = Pi
            Type       = Stationary state frequencies
            Prior      = Dirichlet

      3 --  Parameter  = Alpha
            Type       = Shape of scaled gamma distribution of site rates
            Prior      = Exponential(1.00)

      4 --  Parameter  = Pinvar
            Type       = Proportion of invariable sites
            Prior      = Uniform(0.00,1.00)

      5 --  Parameter  = Ratemultiplier
            Type       = Partition-specific rate multiplier
            Prior      = Fixed(1.0)

      6 --  Parameter  = Tau
            Type       = Topology
            Prior      = All topologies equally probable a priori
            Subparam.  = V

      7 --  Parameter  = V
            Type       = Branch lengths
            Prior      = Unconstrained:GammaDir(1.0,0.1000,1.0,1.0)


   Number of chains per MPI processor = 8

   The MCMC sampler will use the following moves:
      With prob.  Chain will use move
         0.93 %   Dirichlet(Revmat)
         0.93 %   Slider(Revmat)
         0.93 %   Dirichlet(Pi)
         0.93 %   Slider(Pi)
         1.85 %   Multiplier(Alpha)
         1.85 %   Slider(Pinvar)
         9.26 %   ExtSPR(Tau,V)
         9.26 %   ExtTBR(Tau,V)
         9.26 %   NNI(Tau,V)
         9.26 %   ParsSPR(Tau,V)
        37.04 %   Multiplier(V)
        12.96 %   Nodeslider(V)
         5.56 %   TLMultiplier(V)

   Division 1 has 413 unique site patterns
   Initializing conditional likelihoods

   Running benchmarks to automatically select fastest BEAGLE resource... 
   Using BEAGLE v3.1.2 resource 1 for division 1:
      Rsrc Name : Tesla P100-PCIE-16GB
      Impl Name : CUDA-Single
      Flags: PROCESSOR_GPU PRECISION_SINGLE COMPUTATION_SYNCH EIGEN_REAL
             SCALING_MANUAL SCALERS_RAW VECTOR_NONE THREADING_NONE      MODEL STATES: 4
   Initializing invariable-site conditional likelihoods

   Initial log likelihoods and log prior probs for run 1:
      Chain 1 -- -8919.411522 -- 42.620562
      Chain 2 -- -8771.563296 -- 42.620562
      Chain 3 -- -8518.126563 -- 42.620562
      Chain 4 -- -8946.939621 -- 42.620562

   Initial log likelihoods and log prior probs for run 2:
      Chain 1 -- -8592.085448 -- 42.620562
      Chain 2 -- -8440.938886 -- 42.620562
      Chain 3 -- -8688.754041 -- 42.620562
      Chain 4 -- -8894.529587 -- 42.620562


   Using a relative burnin of 25.0 % for diagnostics

   Chain results (20000 generations requested):

      0 -- [-8919.412] (-8771.563) (-8518.127) (-8946.940) * [-8592.085] (-8440.939) (-8688.754) (-8894.530) (...0 remote chains...) 
   1000 -- [-6047.895] (-6056.430) (-6120.001) (-6125.120) * (-6028.749) (-6107.900) [-5983.082] (-6096.173) (...0 remote chains...) -- 0:00:19
   2000 -- (-5966.556) [-5940.409] (-5939.939) (-6091.188) * [-5861.574] (-6038.474) (-5897.188) (-5919.859) (...0 remote chains...) -- 0:00:18
   3000 -- (-5867.629) (-5829.385) [-5814.568] (-5997.532) * (-5844.656) (-5979.460) [-5845.484] (-5852.948) (...0 remote chains...) -- 0:00:17
   4000 -- (-5821.305) [-5746.622] (-5792.056) (-5880.442) * [-5778.875] (-5917.413) (-5817.247) (-5801.371) (...0 remote chains...) -- 0:00:16
   5000 -- (-5782.367) [-5733.577] (-5754.473) (-5811.744) * (-5765.074) (-5849.911) [-5780.278] (-5768.608) (...0 remote chains...) -- 0:00:12

   Average standard deviation of split frequencies: 0.000000

   6000 -- (-5740.663) [-5729.112] (-5736.163) (-5769.044) * [-5731.089] (-5773.169) (-5742.422) (-5735.368) (...0 remote chains...) -- 0:00:11
   7000 -- (-5740.729) (-5735.838) (-5722.699) [-5738.829] * (-5727.399) (-5739.974) [-5733.467] (-5737.091) (...0 remote chains...) -- 0:00:11
   8000 -- (-5732.037) (-5723.254) (-5729.979) [-5733.471] * [-5722.192] (-5733.728) (-5733.738) (-5726.980) (...0 remote chains...) -- 0:00:10
   9000 -- (-5738.576) (-5725.501) (-5730.153) [-5721.744] * [-5725.951] (-5719.629) (-5734.888) (-5722.782) (...0 remote chains...) -- 0:00:09
   10000 -- (-5728.964) (-5739.340) [-5720.749] (-5721.835) * (-5723.985) [-5720.800] (-5727.622) (-5727.411) (...0 remote chains...) -- 0:00:09

   Average standard deviation of split frequencies: 0.001034

   11000 -- [-5729.477] (-5733.644) (-5723.735) (-5721.479) * (-5728.937) (-5724.185) [-5722.981] (-5723.809) (...0 remote chains...) -- 0:00:07
   12000 -- [-5727.612] (-5735.605) (-5720.857) (-5735.472) * (-5733.317) (-5723.410) (-5728.018) [-5723.994] (...0 remote chains...) -- 0:00:06
   13000 -- (-5726.903) [-5720.852] (-5727.963) (-5729.000) * (-5732.313) [-5723.682] (-5723.280) (-5725.235) (...0 remote chains...) -- 0:00:05
   14000 -- (-5734.467) [-5725.775] (-5725.319) (-5739.813) * [-5723.447] (-5734.556) (-5735.317) (-5732.822) (...0 remote chains...) -- 0:00:05
   15000 -- (-5728.179) (-5721.962) (-5724.484) [-5729.233] * [-5725.047] (-5729.352) (-5734.814) (-5729.500) (...0 remote chains...) -- 0:00:04

   Average standard deviation of split frequencies: 0.001378

   16000 -- (-5726.266) (-5719.612) (-5723.972) [-5725.508] * (-5724.286) [-5719.120] (-5726.067) (-5727.467) (...0 remote chains...) -- 0:00:03
   17000 -- (-5724.019) (-5721.611) (-5728.017) [-5732.297] * (-5732.124) [-5721.600] (-5722.171) (-5723.226) (...0 remote chains...) -- 0:00:02
   18000 -- (-5725.945) (-5737.688) [-5727.868] (-5722.582) * [-5727.586] (-5727.064) (-5733.604) (-5724.478) (...0 remote chains...) -- 0:00:01
   19000 -- (-5727.067) (-5721.765) (-5739.332) [-5722.077] * [-5719.814] (-5730.075) (-5734.928) (-5726.458) (...0 remote chains...) -- 0:00:00
   20000 -- (-5732.071) [-5721.275] (-5726.843) (-5717.984) * (-5725.414) (-5733.573) (-5724.125) [-5722.246] (...0 remote chains...) -- 0:00:00

   Average standard deviation of split frequencies: 0.001041

   Continue with analysis? (yes/no): no

   Analysis completed in 17 seconds
   Analysis used 16.29 seconds of CPU time on processor 0
   Likelihood of best state for "cold" chain of run 1 was -5716.00
   Likelihood of best state for "cold" chain of run 2 was -5715.33

   Acceptance rates for the moves in the "cold" chain of run 1:
      With prob.   (last 100)   chain accepted proposals by move
         37.6 %     ( 41 %)     Dirichlet(Revmat)
         58.4 %     ( 57 %)     Slider(Revmat)
         19.3 %     ( 18 %)     Dirichlet(Pi)
         20.3 %     ( 21 %)     Slider(Pi)
         36.3 %     ( 37 %)     Multiplier(Alpha)
         63.5 %     ( 64 %)     Slider(Pinvar)
          0.1 %     (  0 %)     ExtSPR(Tau,V)
          0.2 %     (  0 %)     ExtTBR(Tau,V)
          0.5 %     (  0 %)     NNI(Tau,V)
          0.7 %     (  0 %)     ParsSPR(Tau,V)
         36.5 %     ( 32 %)     Multiplier(V)
         24.2 %     ( 27 %)     Nodeslider(V)
         13.9 %     ( 17 %)     TLMultiplier(V)

   Acceptance rates for the moves in the "cold" chain of run 2:
      With prob.   (last 100)   chain accepted proposals by move
         38.4 %     ( 42 %)     Dirichlet(Revmat)
         54.4 %     ( 56 %)     Slider(Revmat)
         14.6 %     ( 17 %)     Dirichlet(Pi)
         22.1 %     ( 20 %)     Slider(Pi)
         35.1 %     ( 39 %)     Multiplier(Alpha)
         66.7 %     ( 68 %)     Slider(Pinvar)
          0.4 %     (  0 %)     ExtSPR(Tau,V)
          0.3 %     (  0 %)     ExtTBR(Tau,V)
          0.3 %     (  0 %)     NNI(Tau,V)
          0.5 %     (  0 %)     ParsSPR(Tau,V)
         36.6 %     ( 21 %)     Multiplier(V)
         24.9 %     ( 30 %)     Nodeslider(V)
         12.5 %     ( 11 %)     TLMultiplier(V)

   Chain swap information for run 1:

              1     2     3     4 
        --------------------------
      1 |        0.60  0.40  0.26 
      2 |  3285        0.64  0.38 
      3 |  3260  3342        0.58 
      4 |  3403  3348  3362       

   Chain swap information for run 2:

              1     2     3     4 
        --------------------------
      1 |        0.65  0.46  0.30 
      2 |  3451        0.71  0.41 
      3 |  3386  3250        0.57 
      4 |  3304  3271  3338       

   Upper diagonal: Proportion of successful state exchanges between chains
   Lower diagonal: Number of attempted state exchanges between chains

   Chain information:

     ID -- Heat 
    -----------
      1 -- 1.00  (cold chain)
      2 -- 0.91 
      3 -- 0.83 
      4 -- 0.77 

   Heat = 1 / (1 + T * (ID - 1))
      (where T = 0.10 is the temperature and ID is the chain number)


MrBayes > **quit**

   Deleting previously defined characters
   Deleting previously defined taxa
   Quitting program

[user@cn2351 mrbayes]$ **exit**
salloc.exe: Relinquishing job allocation 24838742
[user@biowulf ~]$ 
	
```




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).

 First, create a script, primates-gtr-gamma.nex to execute your MrBayes commands:




```

begin mrbayes;
   set autoclose=yes nowarn=yes;
   execute primates.nex;
   lset nst=6 rates=gamma;
   mcmc nruns=1 ngen=10000 samplefreq=10 file=primates.nex1;
   mcmc file=primates.nex2;
   mcmc file=primates.nex3;
end;

```


 If your run was interrupted, you can resume from the last checkpoint bysetting append=yes with mcmc.



Then create a batch input file (e.g. mrbayes.sh), which differs for CPU and GPU operation.


### CPU


mrbayes.sh:



```

#!/bin/bash
set -e
module load mrbayes
mpirun -np $SLURM_NTASKS mb primates-gtr-gamma.nex

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --no-requeue --ntasks=*#* --ntasks-per-core=1 --constraint x2695 --partition=multinode [--mem-per-cpu=*#*] mrbayes.sh
```


 See [Parallel Jobs in the user guide](/docs/userguide.html#par) for more information.



### GPU


mrbayes.sh:



```

#!/bin/bash
set -e
module load mrbayes
mpirun -np 1 mb primates-gtr-gamma.nex

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --no-requeue --partition=gpu --gres=gpu:*gputype*:1 mrbayes.sh
```


 See [GPU allocation in the user guide](/docs/userguide.html#gpu) for more information.










