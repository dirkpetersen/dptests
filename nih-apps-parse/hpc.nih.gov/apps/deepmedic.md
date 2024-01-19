

document.querySelector('title').textContent = 'funannotate: a pipeline for genome annotation ';
**deepmedic: Deep Learning for segmentation of structures of interest in biomedical 3D scans.** 


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




The deepmedic project aims to offer easy access to Deep Learning 
for segmentation of structures of interest in biomedical 3D scans. 
It is a system that allows the easy creation of a 3D Convolutional Neural Network, 
which can be trained to detect and segment structures 
if corresponding ground truth labels are provided for training. 



### References:


* Konstantinos Kamnitsas, Christian Ledig, Virginia F.J. Newcombe, Joanna P. Simpson, Andrew D. Kane, 
David K. Menon, Daniel Rueckert, Ben Glocker   

*Efficient multi-scale 3D CNN with fully connected CRF for accurate brain lesion segmentation.*   

[Medical Image Analysis
Volume 36, February 2017, Pages 61-78](https://www.sciencedirect.com/science/article/pii/S1361841516301839)


Documentation
* [Deepmedic GitHub page](https://github.com/deepmedic/deepmedic)


Important Notes
* Module Name: deepmedic (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **DM\_HOME**  installation directory
	+ **DM\_BIN**    executables directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
  


```

[user@biowulf]$ **sinteractive --gres=lscratch:10 --mem=10g -c8**
[user@cn0861 ~]$ **module load deepmedic** 
[+] Loading singularity  3.10.5  on cn4185
[+] Loading deepmedic  0.8.4

```

Clone the deepmedic github repository and run deepmedic in CPU mode on the provided example data: 

```

[user@biowulf]$ **wget https://github.com/deepmedic/deepmedic/archive/refs/tags/v0.8.4.tar.gz** 
[user@biowulf]$ **tar -zxf v0.8.4.tar.gz && rm -f v0.8.4.tar.gz && cd deepmedic-0.8.4**
[user@biowulf]$ **deepMedicRun -model ./examples/configFiles/tinyCnn/model/modelConfig.cfg \
 -train examples/configFiles/tinyCnn/train/trainConfigWithValidation.cfg**
Given configuration file:  /vf/users/user/deepmedic/deepmedic-0.8.4/examples/configFiles/tinyCnn/model/modelConfig.cfg
Given configuration file:  /vf/users/user/deepmedic/deepmedic-0.8.4/examples/configFiles/tinyCnn/train/trainConfigWithValidation.cfg
Creating necessary folders for training session...
=============================== logger created =======================================

======================== Starting new session ============================
Command line arguments given:
Namespace(device='cpu', model_cfg='./examples/configFiles/tinyCnn/model/modelConfig.cfg', reset_trainer=False, saved_model=None, test_cfg=None, train_cfg='examples/configFiles/tinyCnn/train/trainConfigWithValidation.cfg')
...
Available devices to Tensorflow:
[name: "/device:CPU:0"
device_type: "CPU"
memory_limit: 268435456
locality {
}
incarnation: 17546255720476806479
, name: "/device:XLA_CPU:0"
device_type: "XLA_CPU"
memory_limit: 17179869184
locality {
}
incarnation: 16465827944861712538
physical_device_desc: "device: XLA_CPU device"
]
CONFIG: The configuration file for the [model] given is: /vf/users/user/deepmedic/deepmedic-0.8.4/examples/configFiles/tinyCnn/model/modelConfig.cfg
=============================================================
========== PARAMETERS FOR MAKING THE ARCHITECTURE ===========
=============================================================
CNN model's name = tinyCnn
~~~~~~~~~~~~~~~~~Model parameters~~~~~~~~~~~~~~~~
Number of Classes (including background) = 5
~Normal Pathway~~
Number of Input Channels = 2
Number of Layers = 3
Number of Feature Maps per layer = [4, 5, 6]
Kernel Dimensions per layer = [[3, 3, 3], [3, 3, 3], [3, 3, 3]]
Padding mode of convs per layer = ['VALID', 'VALID', 'VALID']
Residual connections added at the output of layers (indices from 0) = []
Layers that will be made of Lower Rank (indices from 0) = []
Lower Rank layers will be made of rank = []
~Subsampled Pathway~~
Use subsampled Pathway = True
Number of subsampled pathways that will be built = 1
Number of Layers (per sub-pathway) = [3]
Number of Feature Maps per layer (per sub-pathway) = [[4, 5, 6]]
Kernel Dimensions per layer = [[3, 3, 3], [3, 3, 3], [3, 3, 3]]
Padding mode of convs per layer = ['VALID', 'VALID', 'VALID']
Subsampling Factor per dimension (per sub-pathway) = [[3, 3, 3]]
Residual connections added at the output of layers (indices from 0) = []
Layers that will be made of Lower Rank (indices from 0) = []
Lower Rank layers will be made of rank = []
~Fully Connected Pathway~~
Number of additional FC layers (Excluding the Classif. Layer) = 0
Number of Feature Maps in the additional FC layers = []
Padding mode of convs per layer = ['VALID']
Residual connections added at the output of layers (indices from 0) = []
Layers that will be made of Lower Rank (indices from 0) = []
Dimensions of Kernels in final FC path before classification = [[1, 1, 1]]
~Size Of Image Segments~~
Size of Segments for Training = [25, 25, 25]
Size of Segments for Validation = [7, 7, 7]
Size of Segments for Testing = [45, 45, 45]
~Dropout Rates~~
Drop.R. for each layer in Normal Pathway = []
Drop.R. for each layer in Subsampled Pathway = []
Drop.R. for each layer in FC Pathway (additional FC layers + Classific.Layer at end) = [0.5]
~Weight Initialization~~
Initialization method and params for the conv kernel weights = ['fanIn', 2]
~Activation Function~~
Activation function to use = prelu
~Batch Normalization~~
Apply BN straight on pathways' inputs (eg straight on segments) = [False, False, True]
Batch Normalization uses a rolling average for inference, over this many batches = 60
========== Done with printing session's parameters ==========
=============================================================
CONFIG: The configuration file for the [session] was loaded from: /vf/users/user/deepmedic/deepmedic-0.8.4/examples/configFiles/tinyCnn/train/trainConfigWithValidation.cfg

=============    NEW TRAINING SESSION    ==============


=============================================================
========= PARAMETERS FOR THIS TRAINING SESSION ==============
=============================================================
Session's name = trainSessionWithValidTiny
Model will be loaded from save = None
~Output~~
Main output folder = /vf/users/user/deepmedic/deepmedic-0.8.4/examples/output
Log performance metrics for tensorboard = True
Path and filename to save trained models = /vf/users/user/deepmedic/deepmedic-0.8.4/examples/output/saved_models//trainSessionWithValidTiny//tinyCnn.trainSessionWithValidTiny
~~~~~~~~~~~~~~~~~Generic Information~~~~~~~~~~~~~~~~
Number of Cases for Training = 2
Number of Cases for Validation = 2
~~~~~~~~~~~~~~~~~Training parameters~~~~~~~~~~~~~~~~
Dataframe (csv) filename = None
Filepaths to Channels of the Training Cases = [['/vf/users/user/deepmedic/deepmedic-0.8.4/examples/dataForExamples/brats2015TrainingData/train/brats_2013_pat0005_1/Flair_subtrMeanDivStd.nii.gz', '/vf/users/user/deepmedic/deepmedic-0.8.4/examples/dataForExamples/brats2015TrainingData/train/brats_2013_pat0005_1/T1c_subtrMeanDivStd.nii.gz'], ['/vf/users/user/deepmedic/deepmedic-0.8.4/examples/dataForExamples/brats2015TrainingData/train/brats_2013_pat0006_1/Flair_subtrMeanDivStd.nii.gz', '/vf/users/user/deepmedic/deepmedic-0.8.4/examples/dataForExamples/brats2015TrainingData/train/brats_2013_pat0006_1/T1c_subtrMeanDivStd.nii.gz']]
Filepaths to Ground-Truth labels of the Training Cases = ['/vf/users/user/deepmedic/deepmedic-0.8.4/examples/dataForExamples/brats2015TrainingData/train/brats_2013_pat0005_1/OTMultiClass.nii.gz', '/vf/users/user/deepmedic/deepmedic-0.8.4/examples/dataForExamples/brats2015TrainingData/train/brats_2013_pat0006_1/OTMultiClass.nii.gz']
Filepaths to ROI Masks of the Training Cases = ['/vf/users/user/deepmedic/deepmedic-0.8.4/examples/dataForExamples/brats2015TrainingData/train/brats_2013_pat0005_1/brainmask.nii.gz', '/vf/users/user/deepmedic/deepmedic-0.8.4/examples/dataForExamples/brats2015TrainingData/train/brats_2013_pat0006_1/brainmask.nii.gz']
~ Sampling (train) ~~
Type of Sampling = Fore/Background (0)
Sampling Categories = ['Foreground', 'Background']
Percent of Samples to extract per Sampling Category = [0.5 0.5]
Paths to weight-Maps for sampling of each category = None
~Training Cycle~~
Number of Epochs = 2
Number of Subepochs per epoch = 2
Number of cases to load per Subepoch (for extracting the samples for this subepoch) = 50
Number of Segments loaded per subepoch for Training = 1000. NOTE: This number of segments divided by the batch-size defines the number of optimization-iterations that will be performed every subepoch!
Batch size (train) = 10
Number of parallel processes for sampling = 0
~Learning Rate Schedule~~
Type of schedule = poly
[Predef] Predefined schedule of epochs when the LR will be lowered = None
[Predef] When decreasing Learning Rate, divide LR by = 2.0
[Poly] Initial epochs to wait before lowering LR = 0.6666666666666666
[Poly] Final epoch for the schedule = 2
[Auto] Initial epochs to wait before lowering LR = 5
[Auto] When decreasing Learning Rate, divide LR by = 2.0
[Auto] Minimum increase in validation accuracy (0. to 1.) that resets the waiting counter = 0.0
[Expon] (Deprecated) parameters = {'epochs_wait_before_decr': 0.6666666666666666, 'final_ep_for_sch': 2, 'lr_to_reach_at_last_ep': 0.00390625, 'mom_to_reach_at_last_ep': 0.9}
~Data Augmentation During Training~~
Image level augmentation:
Parameters for image-level augmentation: {'affine': }
 affine: OrderedDict([('prob', 0.7), ('max\_rot\_xyz', (45.0, 45.0, 45.0)), ('max\_scaling', 0.1), ('seed', None), ('interp\_order\_imgs', 1), ('interp\_order\_lbls', 0), ('interp\_order\_roi', 0), ('interp\_order\_wmaps', 1), ('boundary\_mode', 'nearest'), ('cval', 0.0)])
Patch level augmentation:
Mu and std for shift and scale of histograms = {'shift': {'mu': 0.0, 'std': 0.05}, 'scale': {'mu': 1.0, 'std': 0.01}}
Probabilities of reflecting each axis = (0.5, 0.0, 0.0)
Probabilities of rotating planes 0/90/180/270 degrees = {'xy': {'0': 0.8, '90': 0.1, '180': 0.0, '270': 0.1}, 'yz': {'0': 0.0, '90': 0.0, '180': 0.0, '270': 0.0}, 'xz': {'0': 0.0, '90': 0.0, '180': 0.0, '270': 0.0}}
~~~~~~~~~~~~~~~~~Validation parameters~~~~~~~~~~~~~~~~
Perform Validation on Samples throughout training? = True
Perform Full Inference on validation cases every few epochs? = True
Dataframe (csv) filename = None
Filepaths to Channels of Validation Cases = [['/vf/users/user/deepmedic/deepmedic-0.8.4/examples/dataForExamples/brats2015TrainingData/validation/brats\_2013\_pat0003\_1/Flair\_subtrMeanDivStd.nii.gz', '/vf/users/user/deepmedic/deepmedic-0.8.4/examples/dataForExamples/brats2015TrainingData/validation/brats\_2013\_pat0003\_1/T1c\_subtrMeanDivStd.nii.gz'], ['/vf/users/user/deepmedic/deepmedic-0.8.4/examples/dataForExamples/brats2015TrainingData/validation/brats\_2013\_pat0004\_1/Flair\_subtrMeanDivStd.nii.gz', '/vf/users/user/deepmedic/deepmedic-0.8.4/examples/dataForExamples/brats2015TrainingData/validation/brats\_2013\_pat0004\_1/T1c\_subtrMeanDivStd.nii.gz']]
Filepaths to Ground-Truth labels of the Validation Cases = ['/vf/users/user/deepmedic/deepmedic-0.8.4/examples/dataForExamples/brats2015TrainingData/validation/brats\_2013\_pat0003\_1/OTMultiClass.nii.gz', '/vf/users/user/deepmedic/deepmedic-0.8.4/examples/dataForExamples/brats2015TrainingData/validation/brats\_2013\_pat0004\_1/OTMultiClass.nii.gz']
Filepaths to ROI masks for Validation Cases = ['/vf/users/user/deepmedic/deepmedic-0.8.4/examples/dataForExamples/brats2015TrainingData/validation/brats\_2013\_pat0003\_1/brainmask.nii.gz', '/vf/users/user/deepmedic/deepmedic-0.8.4/examples/dataForExamples/brats2015TrainingData/validation/brats\_2013\_pat0004\_1/brainmask.nii.gz']
~~~~~~Validation on Samples throughout Training~~~~~~~
Number of Segments loaded per subepoch for Validation = 5000
Batch size (val on samples) = 50
~ Sampling (val) ~~
Type of Sampling = Uniform (1)
Sampling Categories = ['Uniform']
Percent of Samples to extract per Sampling Category = [1.0]
Paths to weight-maps for sampling of each category = None
~~~~Validation with Full Inference on Validation Cases~~~~~
Perform Full-Inference on Val. cases every that many epochs = 1
Batch size (val on whole volumes) = 10
~Predictions (segmentations and prob maps on val. cases)~~
Save Segmentations = True
Save Probability Maps for each class = [True, True, True, True, True]
Filepaths to save results per case = ['/vf/users/user/deepmedic/deepmedic-0.8.4/examples/output/predictions/trainSessionWithValidTiny/predictions//pred\_brats\_2013\_pat0003\_1.nii.gz', '/vf/users/user/deepmedic/deepmedic-0.8.4/examples/output/predictions/trainSessionWithValidTiny/predictions//pred\_brats\_2013\_pat0004\_1.nii.gz']
Suffixes with which to save segmentations and probability maps = {'segm': 'Segm', 'prob': 'ProbMapClass'}
~Feature Maps~~
Save Feature Maps = False
Min/Max Indices of FMs to visualise per pathway-type and per layer = None
Filepaths to save FMs per case = ['/vf/users/user/deepmedic/deepmedic-0.8.4/examples/output/predictions/trainSessionWithValidTiny/features//pred\_brats\_2013\_pat0003\_1.nii.gz', '/vf/users/user/deepmedic/deepmedic-0.8.4/examples/output/predictions/trainSessionWithValidTiny/features//pred\_brats\_2013\_pat0004\_1.nii.gz']
~Optimization~~
Initial Learning rate = 0.001
Optimizer to use: SGD(0), Adam(1), RmsProp(2) = 2
Parameters for Adam: b1= placeholder, b2=placeholder, e= placeholder
Parameters for RmsProp: rho= 0.9, e= 0.0001
Momentum Type: Classic (0) or Nesterov (1) = 1
Momentum Non-Normalized (0) or Normalized (1) = 1
Momentum Value = 0.6
~Costs~~
Loss functions and their weights = {'xentr': 1.0, 'iou': None, 'dsc': None}
Reweight samples in cost on a per-class basis = {'type': None, 'prms': None, 'schedule': [0, 2]}
L1 Regularization term = 1e-06
L2 Regularization term = 0.0001
~Freeze Weights of Certain Layers~~
Indices of layers from each type of pathway that will be kept fixed (first layer is 0):
Normal pathway's layers to freeze = []
Subsampled pathway's layers to freeze = []
FC pathway's layers to freeze = []
~~~~~~~~~~~~~~~~~ PRE-PROCESSING ~~~~~~~~~~~~~~~~
~Data Compabitibility Checks~~
Check whether input data has correct format (can slow down process) = True
~Padding~~
Pad Input Images = True
~Intensity Normalization~~
Verbosity level = 0
Z-Score parameters = {'apply\_to\_all\_channels': False, 'apply\_per\_channel': None, 'cutoff\_percents': None, 'cutoff\_times\_std': None, 'cutoff\_below\_mean': False}
========== Done with printing session's parameters ==========
=============================================================

=======================================================

=========== Making the CNN graph... ===============
...Building the CNN model...
[Pathway\_NORMAL] is being built...
 Block [0], FMs-In: 2, FMs-Out: 4, Conv Filter dimensions: [3, 3, 3]
WARNING:tensorflow:From /opt/conda/envs/deepmedic/lib/python3.8/site-packages/tensorflow/python/ops/resource\_variable\_ops.py:1659: calling BaseResourceVariable.\_\_init\_\_ (from tensorflow.python.ops.resource\_variable\_ops) with constraint is deprecated and will be removed in a future version.
Instructions for updating:
If using Keras pass \*\_constraint arguments to layers.
 Block [1], FMs-In: 4, FMs-Out: 5, Conv Filter dimensions: [3, 3, 3]
 Block [2], FMs-In: 5, FMs-Out: 6, Conv Filter dimensions: [3, 3, 3]
[Pathway\_SUBSAMPLED[3, 3, 3]] is being built...
 Block [0], FMs-In: 2, FMs-Out: 4, Conv Filter dimensions: [3, 3, 3]
 Block [1], FMs-In: 4, FMs-Out: 5, Conv Filter dimensions: [3, 3, 3]
 Block [2], FMs-In: 5, FMs-Out: 6, Conv Filter dimensions: [3, 3, 3]
[Pathway\_FC] is being built...
 Block [0], FMs-In: 12, FMs-Out: 5, Conv Filter dimensions: [1, 1, 1]
Adding the final Softmax layer...
Finished building the CNN's model.
Pathway [NORMAL], Mode: [train], Input's Shape: (None, 2, 25, 25, 25)
 Block [0], Mode: [train], Input's Shape: (None, 2, 25, 25, 25)
 Block [1], Mode: [train], Input's Shape: (None, 4, 23, 23, 23)
 Block [2], Mode: [train], Input's Shape: (None, 5, 21, 21, 21)
Pathway [NORMAL], Mode: [train], Output's Shape: (None, 6, 19, 19, 19)
Pathway [SUBSAMPLED[3, 3, 3]], Mode: [train], Input's Shape: (None, 2, 13, 13, 13)
 Block [0], Mode: [train], Input's Shape: (None, 2, 13, 13, 13)
 Block [1], Mode: [train], Input's Shape: (None, 4, 11, 11, 11)
 Block [2], Mode: [train], Input's Shape: (None, 5, 9, 9, 9)
Pathway [SUBSAMPLED[3, 3, 3]], Mode: [train], Output's Shape: (None, 6, 7, 7, 7)
Pathway [FC], Mode: [train], Input's Shape: (None, 12, 19, 19, 19)
 Block [0], Mode: [train], Input's Shape: (None, 12, 19, 19, 19)
Pathway [FC], Mode: [train], Output's Shape: (None, 5, 19, 19, 19)
Pathway [NORMAL], Mode: [infer], Input's Shape: (None, 2, 7, 7, 7)
 Block [0], Mode: [infer], Input's Shape: (None, 2, 7, 7, 7)
 Block [1], Mode: [infer], Input's Shape: (None, 4, 5, 5, 5)
 Block [2], Mode: [infer], Input's Shape: (None, 5, 3, 3, 3)
Pathway [NORMAL], Mode: [infer], Output's Shape: (None, 6, 1, 1, 1)
Pathway [SUBSAMPLED[3, 3, 3]], Mode: [infer], Input's Shape: (None, 2, 7, 7, 7)
 Block [0], Mode: [infer], Input's Shape: (None, 2, 7, 7, 7)
 Block [1], Mode: [infer], Input's Shape: (None, 4, 5, 5, 5)
 Block [2], Mode: [infer], Input's Shape: (None, 5, 3, 3, 3)
Pathway [SUBSAMPLED[3, 3, 3]], Mode: [infer], Output's Shape: (None, 6, 1, 1, 1)
Pathway [FC], Mode: [infer], Input's Shape: (None, 12, 1, 1, 1)
 Block [0], Mode: [infer], Input's Shape: (None, 12, 1, 1, 1)
Pathway [FC], Mode: [infer], Output's Shape: (None, 5, 1, 1, 1)
Pathway [NORMAL], Mode: [infer], Input's Shape: (None, 2, 45, 45, 45)
 Block [0], Mode: [infer], Input's Shape: (None, 2, 45, 45, 45)
 Block [1], Mode: [infer], Input's Shape: (None, 4, 43, 43, 43)
 Block [2], Mode: [infer], Input's Shape: (None, 5, 41, 41, 41)
Pathway [NORMAL], Mode: [infer], Output's Shape: (None, 6, 39, 39, 39)
Pathway [SUBSAMPLED[3, 3, 3]], Mode: [infer], Input's Shape: (None, 2, 19, 19, 19)
 Block [0], Mode: [infer], Input's Shape: (None, 2, 19, 19, 19)
 Block [1], Mode: [infer], Input's Shape: (None, 4, 17, 17, 17)
 Block [2], Mode: [infer], Input's Shape: (None, 5, 15, 15, 15)
Pathway [SUBSAMPLED[3, 3, 3]], Mode: [infer], Output's Shape: (None, 6, 13, 13, 13)
Pathway [FC], Mode: [infer], Input's Shape: (None, 12, 39, 39, 39)
 Block [0], Mode: [infer], Input's Shape: (None, 12, 39, 39, 39)
Pathway [FC], Mode: [infer], Output's Shape: (None, 5, 39, 39, 39)
=========== Building Trainer ===========

Building Trainer.
COST: Using cross entropy with weight: 1.0
...Initializing state of the optimizer...
----------- Creating Tensorboard Loggers -----------
Loggers created successfully
-----------=============================-----------
=========== Compiling the Training Function ===========
=======================================================

...Building the training function...
...Collecting ops and feeds for training...
Done.
=========== Compiling the Validation Function =========
...Building the validation function...
...Collecting ops and feeds for validation...
Done.
=========== Compiling the Testing Function ============
...Building the function for testing and visualisation of FMs...
...Collecting ops and feeds for testing...
Done.
=========== Initializing network and trainer variables ===============
All variables were initialized.
Saving the initial model at:/vf/users/user/deepmedic/deepmedic-0.8.4/examples/output/saved\_models//trainSessionWithValidTiny//tinyCnn.trainSessionWithValidTiny.initial.2023-07-26.12.10.28.410496

=======================================================
============== Training the CNN model =================
=======================================================

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~ Starting new Epoch! Epoch #0/2 ~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*
\* Starting new Subepoch: #0/2 \*
\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*
[MAIN|PID:1245871] MULTIPROC: Before Validation in subepoch #0, submitting sampling job for next [VALIDATION].
[VAL|SAMPLER|PID:1245871] :=:=:=:=:=:=: Starting to sample for next [Validation]... :=:=:=:=:=:=:
[VAL|SAMPLER|PID:1245871] Out of [2] subjects given for [Validation], we will sample from maximum [50] per subepoch.
[VAL|SAMPLER|PID:1245871] Shuffled indices of subjects that were randomly chosen: [0, 1]
[VAL|SAMPLER|PID:1245871] Will sample from [2] subjects for next Validation...
[VAL|JOB:0|PID:1245871] Started. (#0/2) sampling job. Load & sample from subject of index (in user's list): 0
[VAL|JOB:0|PID:1245871] Loading subject with 1st channel at: /vf/users/user/deepmedic/deepmedic-0.8.4/examples/dataForExamples/brats2015TrainingData/validation/brats\_2013\_pat0003\_1/Flair\_subtrMeanDivStd.nii.gz
[VAL|JOB:0|PID:1245871] WARN: Loaded labels are dtype [float32]. Rounding and casting to [int16]!
[VAL|JOB:0|PID:1245871] WARN: Loaded ROI-mask is dtype [float64]. Rounding and casting to [int16]!
[VAL|JOB:0|PID:1245871] Done. Samples per category: [Uniform: 2500/2500]
[VAL|JOB:0|PID:1245871] TIMING: [Load: 0.7] [Preproc: 0.1] [Augm-Img: 0.0] [Sample Coords: 0.1] [Extract Sampl: 0.4] [Augm-Samples: 0.0] secs
[VAL|JOB:1|PID:1245871] Started. (#1/2) sampling job. Load & sample from subject of index (in user's list): 1
[VAL|JOB:1|PID:1245871] Loading subject with 1st channel at: /vf/users/user/deepmedic/deepmedic-0.8.4/examples/dataForExamples/brats2015TrainingData/validation/brats\_2013\_pat0004\_1/Flair\_subtrMeanDivStd.nii.gz
[VAL|JOB:1|PID:1245871] WARN: Loaded labels are dtype [float32]. Rounding and casting to [int16]!
[VAL|JOB:1|PID:1245871] WARN: Loaded ROI-mask is dtype [float64]. Rounding and casting to [int16]!
[VAL|JOB:1|PID:1245871] Done. Samples per category: [Uniform: 2500/2500]
[VAL|JOB:1|PID:1245871] TIMING: [Load: 0.6] [Preproc: 0.1] [Augm-Img: 0.0] [Sample Coords: 0.1] [Extract Sampl: 0.4] [Augm-Samples: 0.0] secs
[VAL|SAMPLER|PID:1245871] TIMING: Sampling for next [Validation] lasted: 2.6 secs.
[VAL|SAMPLER|PID:1245871] :=:=:=:=:=:= Finished sampling for next [Validation] =:=:=:=:=:=:
[MAIN|PID:1245871] MULTIPROC: Before Validation in subepoch #0, submitting sampling job for next [TRAINING].
V-V-V-V- Validating for subepoch before starting training iterations -V-V-V-V
[TRA|SAMPLER|PID:1245871] :=:=:=:=:=:=: Starting to sample for next [Training]... :=:=:=:=:=:=:
[VALIDATION] Processed 0/100 batches for this subepoch...
[TRA|SAMPLER|PID:1245871] Out of [2] subjects given for [Training], we will sample from maximum [50] per subepoch.
[TRA|SAMPLER|PID:1245871] Shuffled indices of subjects that were randomly chosen: [0, 1]
[TRA|SAMPLER|PID:1245871] Will sample from [2] subjects for next Training...
[TRA|JOB:0|PID:1245871] Started. (#0/2) sampling job. Load & sample from subject of index (in user's list): 0
[TRA|JOB:0|PID:1245871] Loading subject with 1st channel at: /vf/users/user/deepmedic/deepmedic-0.8.4/examples/dataForExamples/brats2015TrainingData/train/brats\_2013\_pat0005\_1/Flair\_subtrMeanDivStd.nii.gz
[VALIDATION] Processed 20/100 batches for this subepoch...
[TRA|JOB:0|PID:1245871] WARN: Loaded labels are dtype [float32]. Rounding and casting to [int16]!
[VALIDATION] Processed 40/100 batches for this subepoch...
[VALIDATION] Processed 60/100 batches for this subepoch...
[TRA|JOB:0|PID:1245871] WARN: Loaded ROI-mask is dtype [float64]. Rounding and casting to [int16]!
[VALIDATION] Processed 80/100 batches for this subepoch...
[VALIDATION] Processed 100/100 batches for this subepoch...

+++++++++++++++++++++++ Reporting Accuracy over whole subepoch +++++++++++++++++++++++
VALIDATION: Epoch #0, Subepoch #0, Overall: mean accuracy: 0.0436 => Correctly-Classified-Voxels/All-Predicted-Voxels = 218/5000
+++++++++++++++ Reporting Accuracy over whole subepoch for Class-0 ++++++++ [Whole Foreground (Pos) Vs Background (Neg)] ++++++++++++++++
VALIDATION: Epoch #0, Subepoch #0, Class-0: mean accuracy: 0.0950 => (TruePos+TrueNeg)/All-Predicted-Voxels = 475/5000
VALIDATION: Epoch #0, Subepoch #0, Class-0: mean sensitivity: 1.0000 => TruePos/RealPos = 475/475
...
VALIDATION: Epoch #0, Subepoch #0, Class-4: mean Dice: 0.0258
=============== LOGGING TO TENSORBOARD ===============
Logging VALIDATION metrics
Epoch: 0 | Subepoch 0
Step number (index of subepoch since start): 0
--- Logging per class metrics ---
Logged metrics: ['samples: accuracy', 'samples: sensitivity', 'samples: precision', 'samples: specificity', 'samples: Dice']
======================================================
TIMING: Validation on batches of subepoch #0 lasted: 0.9 secs.
[TRA|JOB:0|PID:1245871] Done. Samples per category: [Foreground: 250/250] [Background: 250/250]
[TRA|JOB:0|PID:1245871] TIMING: [Load: 0.8] [Preproc: 0.1] [Augm-Img: 3.9] [Sample Coords: 0.1] [Extract Sampl: 0.1] [Augm-Samples: 0.3] secs
[TRA|JOB:1|PID:1245871] Started. (#1/2) sampling job. Load & sample from subject of index (in user's list): 1
[TRA|JOB:1|PID:1245871] Loading subject with 1st channel at: /vf/users/user/deepmedic/deepmedic-0.8.4/examples/dataForExamples/brats2015TrainingData/train/brats\_2013\_pat0006\_1/Flair\_subtrMeanDivStd.nii.gz
[TRA|JOB:1|PID:1245871] WARN: Loaded labels are dtype [float32]. Rounding and casting to [int16]!
[TRA|JOB:1|PID:1245871] WARN: Loaded ROI-mask is dtype [float64]. Rounding and casting to [int16]!
[TRA|JOB:1|PID:1245871] Done. Samples per category: [Foreground: 250/250] [Background: 250/250]
[TRA|JOB:1|PID:1245871] TIMING: [Load: 0.6] [Preproc: 0.1] [Augm-Img: 0.0] [Sample Coords: 0.1] [Extract Sampl: 0.1] [Augm-Samples: 0.3] secs
[TRA|SAMPLER|PID:1245871] TIMING: Sampling for next [Training] lasted: 6.8 secs.
[TRA|SAMPLER|PID:1245871] :=:=:=:=:=:= Finished sampling for next [Training] =:=:=:=:=:=:
[MAIN|PID:1245871] MULTIPROC: Before Training in subepoch #0, submitting sampling job for next [VALIDATION].
-T-T-T-T- Training for this subepoch... May take a few minutes... -T-T-T-T-
[VAL|SAMPLER|PID:1245871] :=:=:=:=:=:=: Starting to sample for next [Validation]... :=:=:=:=:=:=:
[TRAINING] Processed 0/100 batches for this subepoch...
[VAL|SAMPLER|PID:1245871] Out of [2] subjects given for [Validation], we will sample from maximum [50] per subepoch.
[VAL|SAMPLER|PID:1245871] Shuffled indices of subjects that were randomly chosen: [0, 1]
[VAL|SAMPLER|PID:1245871] Will sample from [2] subjects for next Validation...
[VAL|JOB:0|PID:1245871] Started. (#0/2) sampling job. Load & sample from subject of index (in user's list): 0
[VAL|JOB:0|PID:1245871] Loading subject with 1st channel at: /vf/users/user/deepmedic/deepmedic-0.8.4/examples/dataForExamples/brats2015TrainingData/validation/brats\_2013\_pat0003\_1/Flair\_subtrMeanDivStd.nii.gz
[VAL|JOB:0|PID:1245871] WARN: Loaded labels are dtype [float32]. Rounding and casting to [int16]!
[VAL|JOB:0|PID:1245871] WARN: Loaded ROI-mask is dtype [float64]. Rounding and casting to [int16]!
[VAL|JOB:0|PID:1245871] Done. Samples per category: [Uniform: 2500/2500]
[VAL|JOB:0|PID:1245871] TIMING: [Load: 0.6] [Preproc: 0.1] [Augm-Img: 0.0] [Sample Coords: 0.0] [Extract Sampl: 0.5] [Augm-Samples: 0.0] secs
[VAL|JOB:1|PID:1245871] Started. (#1/2) sampling job. Load & sample from subject of index (in user's list): 1
[VAL|JOB:1|PID:1245871] Loading subject with 1st channel at: /vf/users/user/deepmedic/deepmedic-0.8.4/examples/dataForExamples/brats2015TrainingData/validation/brats\_2013\_pat0004\_1/Flair\_subtrMeanDivStd.nii.gz
[VAL|JOB:1|PID:1245871] WARN: Loaded labels are dtype [float32]. Rounding and casting to [int16]!
[VAL|JOB:1|PID:1245871] WARN: Loaded ROI-mask is dtype [float64]. Rounding and casting to [int16]!
[VAL|JOB:1|PID:1245871] Done. Samples per category: [Uniform: 2500/2500]
[VAL|JOB:1|PID:1245871] TIMING: [Load: 0.7] [Preproc: 0.2] [Augm-Img: 0.0] [Sample Coords: 0.1] [Extract Sampl: 0.6] [Augm-Samples: 0.0] secs
[VAL|SAMPLER|PID:1245871] TIMING: Sampling for next [Validation] lasted: 3.0 secs.
[VAL|SAMPLER|PID:1245871] :=:=:=:=:=:= Finished sampling for next [Validation] =:=:=:=:=:=:
[TRAINING] Processed 20/100 batches for this subepoch...
[TRAINING] Processed 40/100 batches for this subepoch...
[TRAINING] Processed 60/100 batches for this subepoch...
[TRAINING] Processed 80/100 batches for this subepoch...
[TRAINING] Processed 100/100 batches for this subepoch...

+++++++++++++++++++++++ Reporting Accuracy over whole subepoch +++++++++++++++++++++++
TRAINING: Epoch #0, Subepoch #0, Overall: mean accuracy: 0.3992 => Correctly-Classified-Voxels/All-Predicted-Voxels = 2738260/6859000
TRAINING: Epoch #0, Subepoch #0, Overall: mean cost: 1.38944
+++++++++++++++ Reporting Accuracy over whole subepoch for Class-0 ++++++++ [Whole Foreground (Pos) Vs Background (Neg)] ++++++++++++++++
...
TRAINING: Epoch #0, Subepoch #0, Class-4: mean precision: 0.3182 => TruePos/(TruePos+FalsePos) = 299811/942154
TRAINING: Epoch #0, Subepoch #0, Class-4: mean specificity: 0.8967 => TrueNeg/RealNeg = 5577735/6220078
TRAINING: Epoch #0, Subepoch #0, Class-4: mean Dice: 0.3792
=============== LOGGING TO TENSORBOARD ===============
Logging TRAINING metrics
Epoch: 0 | Subepoch 0
Step number (index of subepoch since start): 0
--- Logging average metrics for all classes ---
Logged metrics: ['samples: accuracy', 'samples: cost']
--- Logging per class metrics ---
Logged metrics: ['samples: accuracy', 'samples: sensitivity', 'samples: precision', 'samples: specificity', 'samples: Dice']
======================================================
TIMING: Training on batches of this subepoch #0 lasted: 16.7 secs.

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*
\* Starting new Subepoch: #1/2 \*
\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*
[MAIN|PID:1245871] MULTIPROC: Before Validation in subepoch #1, submitting sampling job for next [TRAINING].
V-V-V-V- Validating for subepoch before starting training iterations -V-V-V-V
[TRA|SAMPLER|PID:1245871] :=:=:=:=:=:=: Starting to sample for next [Training]... :=:=:=:=:=:=:
[VALIDATION] Processed 0/100 batches for this subepoch...
[TRA|SAMPLER|PID:1245871] Out of [2] subjects given for [Training], we will sample from maximum [50] per subepoch.
[TRA|SAMPLER|PID:1245871] Shuffled indices of subjects that were randomly chosen: [0, 1]
[TRA|SAMPLER|PID:1245871] Will sample from [2] subjects for next Training...
[TRA|JOB:0|PID:1245871] Started. (#0/2) sampling job. Load & sample from subject of index (in user's list): 0
[TRA|JOB:0|PID:1245871] Loading subject with 1st channel at: /vf/users/user/deepmedic/deepmedic-0.8.4/examples/dataForExamples/brats2015TrainingData/train/brats\_2013\_pat0005\_1/Flair\_subtrMeanDivStd.nii.gz
[VALIDATION] Processed 20/100 batches for this subepoch...
[VALIDATION] Processed 40/100 batches for this subepoch...
[VALIDATION] Processed 60/100 batches for this subepoch...
[VALIDATION] Processed 80/100 batches for this subepoch...
[VALIDATION] Processed 100/100 batches for this subepoch...

+++++++++++++++++++++++ Reporting Accuracy over whole subepoch +++++++++++++++++++++++
VALIDATION: Epoch #0, Subepoch #1, Overall: mean accuracy: 0.8306 => Correctly-Classified-Voxels/All-Predicted-Voxels = 4153/5000
...
VALIDATION: Epoch #0, Subepoch #1, Class-4: mean sensitivity: 0.7910 => TruePos/RealPos = 53/67
VALIDATION: Epoch #0, Subepoch #1, Class-4: mean precision: 0.1312 => TruePos/(TruePos+FalsePos) = 53/404
VALIDATION: Epoch #0, Subepoch #1, Class-4: mean specificity: 0.9288 => TrueNeg/RealNeg = 4582/4933
VALIDATION: Epoch #0, Subepoch #1, Class-4: mean Dice: 0.2251
=============== LOGGING TO TENSORBOARD ===============
Logging VALIDATION metrics
Epoch: 0 | Subepoch 1
Step number (index of subepoch since start): 1
--- Logging per class metrics ---
Logged metrics: ['samples: accuracy', 'samples: sensitivity', 'samples: precision', 'samples: specificity', 'samples: Dice']
======================================================
TIMING: Validation on batches of subepoch #1 lasted: 0.6 secs.
[TRA|JOB:0|PID:1245871] WARN: Loaded labels are dtype [float32]. Rounding and casting to [int16]!
[TRA|JOB:0|PID:1245871] WARN: Loaded ROI-mask is dtype [float64]. Rounding and casting to [int16]!
[TRA|JOB:0|PID:1245871] Done. Samples per category: [Foreground: 250/250] [Background: 250/250]
...
TIMING: Training process lasted: 134.6 secs.
Closing worker pool.
 Saving the final model at:/vf/users/userXC
 /deepmedic/deepmedic-0.8.4/examples/output/saved\_models//trainSessionWithValidTiny//tinyCnn.trainSessionWithValidTiny.final.2023-07-26.12.12.43.174418
The whole do\_training() function has finished.

=======================================================
=========== Training session finished =================
=======================================================
Finished.

```

In order to run deepmedic in GPU mode on the same data, pass it an additional option -dev cuda0:

```

[user@biowulf]$ **deepMedicRun -model ./examples/configFiles/tinyCnn/model/modelConfig.cfg \
 -train examples/configFiles/tinyCnn/train/trainConfigWithValidation.cfgi -dev cuda0** 
...
TIMING: Training process lasted: 85.7 secs.
Closing worker pool.
Saving the final model at:/vf/users/$USER/deepmedic/deepmedic-0.8.4/examples/output/saved_models//trainSessionWithValidTiny//tinyCnn.trainSessionWithValidTiny.final.2023-07-26.12.19.06.158002
The whole do_training() function has finished.

=======================================================
=========== Training session finished =================
=======================================================
Finished.
[user@cn0861 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226

```





