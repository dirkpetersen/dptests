

document.querySelector('title').textContent = 'cellpose on Biowulf';
cellpose on Biowulf


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



The cellpose is a generalist algorithm for cellular segmentation with human-in-the-loop capabilities.






### References:


* Pachitariu, M. & Stringer, C. *Cellpose 2.0: how to train your own model.* Nature methods.2022

 [Journal](https://www.nature.com/articles/s41592-022-01663-4)


Documentation
* cellpose Github:[Github](https://github.com/MouseLand/cellpose)


Important Notes
* Module Name: cellpose (see [the modules page](/apps/modules.html) for more information)
 * The cellpose command lines could be run as:
 
```

	cellpose -h
	
```
* cellpose GUI could to be opened through [NoMachine](https://hpc.nih.gov/docs/nx.html): 
 
```

       cellpose-gui
       
```
* cellpose GUI could also be opened through visual partition [svis](https://hpc.nih.gov/docs/svis.html) desktop session:
 
```

       vglrun cellpose-gui
       
```


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --gres=gpu:p100:1 --mem=8g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load cellpose**
[user@cn3144 ~]$ **mkdir /data/$USER/cellpose\_test/**
[user@cn3144 ~]$ **cd /data/$USER/cellpose\_test/**
[user@cn3144 ~]$ **cellpose -h**
/opt/conda/lib/python3.8/site-packages/pyqtgraph/colors/palette.py:1: RuntimeWarning: PyQtGraph supports Qt version >= 5.15, but 5.12.9 detected.
  from ..Qt import QtGui
usage: cellpose [-h] [--use_gpu] [--gpu_device GPU_DEVICE] [--check_mkl] [--dir DIR] [--image_path IMAGE_PATH] [--look_one_level_down]
                [--img_filter IMG_FILTER] [--channel_axis CHANNEL_AXIS] [--z_axis Z_AXIS] [--chan CHAN] [--chan2 CHAN2] [--invert] [--all_channels]
                [--pretrained_model PRETRAINED_MODEL] [--add_model ADD_MODEL] [--unet] [--nclasses NCLASSES] [--no_resample] [--net_avg] [--no_interp]
                [--no_norm] [--do_3D] [--diameter DIAMETER] [--stitch_threshold STITCH_THRESHOLD] [--fast_mode] [--flow_threshold FLOW_THRESHOLD]
                [--cellprob_threshold CELLPROB_THRESHOLD] [--anisotropy ANISOTROPY] [--exclude_on_edges] [--save_png] [--save_tif] [--no_npy]
                [--savedir SAVEDIR] [--dir_above] [--in_folders] [--save_flows] [--save_outlines] [--save_ncolor] [--save_txt] [--train] [--train_size]
                [--test_dir TEST_DIR] [--mask_filter MASK_FILTER] [--diam_mean DIAM_MEAN] [--learning_rate LEARNING_RATE] [--weight_decay WEIGHT_DECAY]
                [--n_epochs N_EPOCHS] [--batch_size BATCH_SIZE] [--min_train_masks MIN_TRAIN_MASKS] [--residual_on RESIDUAL_ON] [--style_on STYLE_ON]
                [--concatenation CONCATENATION] [--save_every SAVE_EVERY] [--save_each] [--verbose]

cellpose parameters

optional arguments:
  -h, --help            show this help message and exit
  --verbose             show information about running and settings and save to log

hardware arguments:
  --use_gpu             use gpu if torch with cuda installed
  --gpu_device GPU_DEVICE
                        which gpu device to use, use an integer for torch, or mps for M1
  --check_mkl           check if mkl working

input image arguments:
  --dir DIR             folder containing data to run or train on.
  --image_path IMAGE_PATH
                        if given and --dir not given, run on single image instead of folder (cannot train with this option)
  --look_one_level_down
                        run processing on all subdirectories of current folder
  --img_filter IMG_FILTER
                        end string for images to run on
  --channel_axis CHANNEL_AXIS
                        axis of image which corresponds to image channels
  --z_axis Z_AXIS       axis of image which corresponds to Z dimension
  --chan CHAN           channel to segment; 0: GRAY, 1: RED, 2: GREEN, 3: BLUE. Default: 0
  --chan2 CHAN2         nuclear channel (if cyto, optional); 0: NONE, 1: RED, 2: GREEN, 3: BLUE. Default: 0
  --invert              invert grayscale channel
  --all_channels        use all channels in image if using own model and images with special channels

model arguments:
  --pretrained_model PRETRAINED_MODEL
                        model to use for running or starting training
  --add_model ADD_MODEL
                        model path to copy model to hidden .cellpose folder for using in GUI/CLI
  --unet                run standard unet instead of cellpose flow output
  --nclasses NCLASSES   if running unet, choose 2 or 3; cellpose always uses 3

algorithm arguments:
  --no_resample         disable dynamics on full image (makes algorithm faster for images with large diameters)
  --net_avg             run 4 networks instead of 1 and average results
  --no_interp           do not interpolate when running dynamics (was default)
  --no_norm             do not normalize images (normalize=False)
  --do_3D               process images as 3D stacks of images (nplanes x nchan x Ly x Lx
  --diameter DIAMETER   cell diameter, if 0 will use the diameter of the training labels used in the model, or with built-in model will estimate diameter
                        for each image
  --stitch_threshold STITCH_THRESHOLD
                        compute masks in 2D then stitch together masks with IoU>0.9 across planes
  --fast_mode           now equivalent to --no_resample; make code run faster by turning off resampling
  --flow_threshold FLOW_THRESHOLD
                        flow error threshold, 0 turns off this optional QC step. Default: 0.4
  --cellprob_threshold CELLPROB_THRESHOLD
                        cellprob threshold, default is 0, decrease to find more and larger masks
  --anisotropy ANISOTROPY
                        anisotropy of volume in 3D
  --exclude_on_edges    discard masks which touch edges of image

output arguments:
  --save_png            save masks as png and outlines as text file for ImageJ
  --save_tif            save masks as tif and outlines as text file for ImageJ
  --no_npy              suppress saving of npy
  --savedir SAVEDIR     folder to which segmentation results will be saved (defaults to input image directory)
  --dir_above           save output folders adjacent to image folder instead of inside it (off by default)
  --in_folders          flag to save output in folders (off by default)
  --save_flows          whether or not to save RGB images of flows when masks are saved (disabled by default)
  --save_outlines       whether or not to save RGB outline images when masks are saved (disabled by default)
  --save_ncolor         whether or not to save minimal "n-color" masks (disabled by default
  --save_txt            flag to enable txt outlines for ImageJ (disabled by default)

training arguments:
  --train               train network using images in dir
  --train_size          train size network at end of training
  --test_dir TEST_DIR   folder containing test data (optional)
  --mask_filter MASK_FILTER
                        end string for masks to run on. use "_seg.npy" for manual annotations from the GUI. Default: _masks
  --diam_mean DIAM_MEAN
                        mean diameter to resize cells to during training -- if starting from pretrained models it cannot be changed from 30.0
  --learning_rate LEARNING_RATE
                        learning rate. Default: 0.2
  --weight_decay WEIGHT_DECAY
                        weight decay. Default: 1e-05
  --n_epochs N_EPOCHS   number of epochs. Default: 500
  --batch_size BATCH_SIZE
                        batch size. Default: 8
  --min_train_masks MIN_TRAIN_MASKS
                        minimum number of masks a training image must have to be used. Default: 5
  --residual_on RESIDUAL_ON
                        use residual connections
  --style_on STYLE_ON   use style vector
  --concatenation CONCATENATION
                        concatenate downsampled layers with upsampled layers (off by default which means they are added)
  --save_every SAVE_EVERY
                        number of epochs to skip between saves. Default: 100
  --save_each           save the model under a different filename per --save_every epoch for later comparsion


[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


 Open cellpose-gui inside of NoMachine:

```

[user@biowulf]$ **sinteractive --gres=gpu:p100:1 --mem=8g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load cellpose**
[user@cn3144 ~]$ **cellpose-gui**

```


![cellpose-gui window](/images/cellpose.png)

 Open cellpose-gui through visual partition [svis](https://hpc.nih.gov/docs/svis.html) desktop session:

```

[user@cn0655]$ **module load virtualgl**
[user@cn0655]$ **module load cellpose**
[user@cn0655]$ **vglrun cellpose-gui**

```


![cellpose-gui window](/images/cellpose.png)


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. cellpose.sh). For example:




hljs.highlightAll();

```


#!/bin/bash
set -e
module load cellpose
cellpose --train --dir /path/to/images/ --pretrained_model None --chan 3 --chan2 0 --use_gpu --verbose


```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --partition=gpu --gres=gpu:p100:1 --mem=8g cellpose.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. cellpose.swarm). For example:



```


cellpose --train --dir /path/to/images1/ --pretrained_model None --chan 3 --chan2 0 --use_gpu --verbose 
cellpose --train --dir /path/to/images2/ --pretrained_model None --chan 3 --chan2 0 --use_gpu --verbose

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f cellpose.swarm [-t #] [-g #] --partition=gpu --gres=gpu:p100:1 --module cellpose
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module cellpose Loads the cellpose module for each subjob in the swarm
 | |
 | |
 | |














