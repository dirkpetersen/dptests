

document.querySelector('title').textContent = 'elastix on Biowulf ';
elastix on Biowulf 



|  |
| --- |
| 
Quick Links
[Interactive job on Biowulf](#int)
[Batch job on Biowulf](#serial)
[Swarm of jobs](#swarm)
 |



Description

elastix is a tool for (medical) image registration. It's based on the Insight
Segmentation and Registration Toolkit (ITK).



There may be multiple versions of elastix available. An easy way of selecting the
version is to use [modules](/apps/modules.html). To see the modules
available, type



```

module avail elastix 

```

To select a module use



```

module load elastix/[version]

```

where `[version]` is the version of choice.


### Environment variables set


* `$PATH`


### References


* S. Klein, M. Staring, K. Murphy, M.A. Viergever, J.P.W. Pluim. 
 *elastix: a toolbox for intensity based medical image registration*. 
 IEEE Transactions on Medical Imaging 2010, 29:196-205.
 [Pubmed](https://www.ncbi.nlm.nih.gov/pubmed/19923044) | 
 PMC  | 
 [Journal](http://ieeexplore.ieee.org/document/5338015/)


### Documentation


* [Home page](http://elastix.isi.uu.nl/index.php)
* [Github page](https://github.com/SuperElastix/elastix)




Interactive job on Biowulf
Allocate an interactive session with [sinteractive](/docs/userguide.html#int)
and use as shown below



```

biowulf$ **sinteractive** 
node$ **module load elastix**
[+] Loading elastix 4.9
node$ **elastix --help**
elastix version: 4.900

elastix registers a moving image to a fixed image.
The registration-process is specified in the parameter file.
  --help, -h displays this message and exit
  --version  output version information and exit

Call elastix from the command line with mandatory arguments:
  -f        fixed image
  -m        moving image
  -out      output directory
  -p        parameter file, elastix handles 1 or more "-p"

Optional extra commands:
  -fMask    mask for fixed image
  -mMask    mask for moving image
  -t0       parameter file for initial transform
  -priority set the process priority to high, abovenormal, normal (default),
            belownormal, or idle (Windows only option)
  -threads  set the maximum number of threads of elastix

The parameter-file must contain all the information necessary for elastix to run properly. That includes which metric to use, which optimizer, which transform, etc. It must also contain information specific for the metric, optimizer, transform, etc. For a usable parameter-file, see the website.

Need further help?
Check the website http://elastix.isi.uu.nl, or mail elastix@bigr.nl.
node$ **mkdir res**
node$ **cp $ELASTIX\_TEST/\* .**
node$ **elastix -f tile\_A.tif -m tile\_B.tif -out res -p parameters.2D.txt** 
...
Running elastix with parameter file 0: "parameters.2D.txt".

Current time: Tue May 15 09:10:28 2018.
Reading the elastix parameters from file ...
...
Resolution: 0
...
Elastix initialization of all components (for this resolution) took: 9 ms.
Initialization of TransformBendingEnergy metric took: 0 ms.
Starting automatic parameter estimation for AdaptiveStochasticGradientDescent ...
WARNING: The parameter "ASGDParameterEstimationMethod", requested at entry number 0, does not exist at all.
  The default value "Original" is used instead.
  Computing JacobianTerms ...
  Computing the Jacobian terms took 0.083924s
  NumberOfGradientMeasurements to estimate sigma_i: 2
  Sampling gradients ...
  Sampling the gradients took 0.258976s
Automatic parameter estimation took 0.35s
1:ItNr  2:Metric        3a:Time 3b:StepSize     4:||Gradient||  Time[ms]
0       0.000000        0.000000        0.000000        0.000000        373.1
1       0.000000        0.000000        0.000000        0.000000        5.5
2       0.000000        0.000000        0.000000        0.000000        5.6
3       0.000000        0.000000        0.000000        0.000000        5.6
4       0.000000        0.000000        0.000000        0.000000        5.5
...
96      0.000000        0.000000        0.000000        0.000000        5.4
97      0.000000        0.000000        0.000000        0.000000        5.4
98      0.000000        0.000000        0.000000        0.000000        5.3
99      0.000000        0.000000        0.000000        0.000000        5.5
Time spent in resolution 0 (ITK initialization and iterating): 0.917 s.
Stopping condition: Maximum number of iterations has been reached.
Settings of AdaptiveStochasticGradientDescent in resolution 0:
( SP_a 0.000000 )
( SP_A 20.000000 )
( SP_alpha 1.000000 )
( SigmoidMax 1.000000 )
( SigmoidMin -0.990000 )
( SigmoidScale 0.000000 )


Resolution: 1
...
Elastix initialization of all components (for this resolution) took: 20 ms.
Initialization of TransformBendingEnergy metric took: 0 ms.
Starting automatic parameter estimation for AdaptiveStochasticGradientDescent ...
WARNING: The parameter "ASGDParameterEstimationMethod", requested at entry number 0, does not exist at all.
  The default value "Original" is used instead.
  Computing JacobianTerms ...
  Computing the Jacobian terms took 0.629961s
  NumberOfGradientMeasurements to estimate sigma_i: 2
  Sampling gradients ...
  Sampling the gradients took 0.246773s
Automatic parameter estimation took 0.89s
1:ItNr  2:Metric        3a:Time 3b:StepSize     4:||Gradient||  Time[ms]
0       0.000000        0.000000        50198347.615187 0.000000        908.2
1       0.000000        0.000000        50198347.615187 0.000000        7.1
2       0.000000        0.000000        50198347.615187 0.000000        7.1
3       0.000000        0.000000        50198347.615187 0.000000        7.0
...
95      0.000000        0.000000        50198347.615187 0.000000        6.9
96      0.000000        0.000000        50198347.615187 0.000000        7.0
97      0.000000        0.000000        50198347.615187 0.000000        6.9
98      0.000000        0.000000        50198347.615187 0.000000        7.0
99      0.000000        0.000000        50198347.615187 0.000000        7.0
Time spent in resolution 1 (ITK initialization and iterating): 1.606 s.
Stopping condition: Maximum number of iterations has been reached.
Settings of AdaptiveStochasticGradientDescent in resolution 1:
( SP_a 1054165299.918932 )
( SP_A 20.000000 )
( SP_alpha 1.000000 )
( SigmoidMax 1.000000 )
( SigmoidMin -0.722400 )
( SigmoidScale 0.000000 )


Resolution: 2
...
Elastix initialization of all components (for this resolution) took: 54 ms.
Initialization of TransformBendingEnergy metric took: 5 ms.
Starting automatic parameter estimation for AdaptiveStochasticGradientDescent ...
WARNING: The parameter "ASGDParameterEstimationMethod", requested at entry number 0, does not exist at all.
  The default value "Original" is used instead.
  Computing JacobianTerms ...
  Computing the Jacobian terms took 9.428071s
  NumberOfGradientMeasurements to estimate sigma_i: 2
  Sampling gradients ...
  Sampling the gradients took 0.347363s
Automatic parameter estimation took 9.79s
1:ItNr  2:Metric        3a:Time 3b:StepSize     4:||Gradient||  Time[ms]
0       0.000000        0.000000        4335379.536497  0.000000        9830.5
1       0.000000        0.000000        4335379.536497  0.000000        21.0
2       0.000000        0.000000        4335379.536497  0.000000        21.0
3       0.000000        0.000000        4335379.536497  0.000000        20.8
...
96      0.000000        0.000000        4335379.536497  0.000000        20.9
97      0.000000        0.000000        4335379.536497  0.000000        20.8
98      0.000000        0.000000        4335379.536497  0.000000        21.0
99      0.000000        0.000000        4335379.536497  0.000000        21.0
Time spent in resolution 2 (ITK initialization and iterating): 11.913 s.
Stopping condition: Maximum number of iterations has been reached.
Settings of AdaptiveStochasticGradientDescent in resolution 2:
( SP_a 91042970.266438 )
( SP_A 20.000000 )
( SP_alpha 1.000000 )
( SigmoidMax 1.000000 )
( SigmoidMin -0.902823 )
( SigmoidScale 0.000000 )



Creating the TransformParameterFile took 0.15s
Registration result checksum: 2802733262

Skipping applying final transform, no resulting output image generated.

Final metric value  = 0.000000
Settings of AdaptiveStochasticGradientDescent for all resolutions:
( SP_a 0.000000 1054165299.918932 91042970.266438 )
( SP_A 20.000000 20.000000 20.000000 )
( SP_alpha 1.000000 1.000000 1.000000 )
( SigmoidMax 1.000000 1.000000 1.000000 )
( SigmoidMin -0.990000 -0.722400 -0.902823 )
( SigmoidScale 0.000000 0.000000 0.000000 )

Time spent on saving the results, applying the final transform etc.: 161 ms.
Running elastix with parameter file 0: "parameters.2D.txt", has finished.


Current time: Tue May 15 09:10:44 2018.
Time used for running elastix with this parameter file:
  15.8s.

-------------------------------------------------------------------------

Total time elapsed: 15.9s.

node$ **exit**
biowulf$

```



Batch job on Biowulf
Create a batch script similar to the following example:



```

#! /bin/bash
# this file is elastix.sh

module load elastix || exit 1
elastix -f fixedImage.ext \
        -m movingImage.ext \
        -out outputDirectory \
        -p parameterFile.txt

```

Submit to the queue with [sbatch](/docs/userguide.html):



```

biowulf$ **sbatch elastix.batch**

```



Swarm of jobs on Biowulf
Create a swarm command file similar to the following example:



```

# this file is elastix.swarm
elastix -f fixedImage.ext -m movingImage.ext -out outputDirectory1 -p parameterFile1.txt
elastix -f fixedImage.ext -m movingImage.ext -out outputDirectory2 -p parameterFile2.txt
elastix -f fixedImage.ext -m movingImage.ext -out outputDirectory3 -p parameterFile3.txt

```

And submit to the queue with [swarm](/apps/swarm.html)



```

biowulf$ **swarm -f elastix.swarm**

```



