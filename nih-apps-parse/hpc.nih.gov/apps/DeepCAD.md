

document.querySelector('title').textContent = 'DeepCAD: Deep self-supervised learning for calcium imaging denoising..';
**DeepCAD: Deep self-supervised learning for calcium imaging denoising..**


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



DeepCAD is a self-supervised deep-learning method 
for spatiotemporal enhancement of calcium imaging data 
that does not require any high signal-to-noise ratio (SNR) observations. 
DeepCAD suppresses detection noise and improves the SNR more than tenfold, 
which reinforces the accuracy of neuron extraction and spike inference 
and facilitates the functional analysis of neural circuits.



### References:


 Xinyang Li, Guoxun Zhang, Jiamin Wu, Yuanlong Zhang, Zhifeng Zhao,
Xing Lin, Hui Qiao, Hao Xie, Haoqian Wang, Lu Fang and Qionghai Dai   

*Reinforcing neuron extraction and spike inference
in calcium imaging using deep self-supervised denoising*    

[Nature Methods (2021)](https://www.nature.com/articles/s41592-021-01225-0), VOL 18, pp. 1395–1400.


Documentation
* [DeepCAD Github page](https://github.com/cabooster/DeepCAD)


Important Notes
* Module Name: DeepCAD (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **DEEPCAD\_HOME**  installation directory
	+ **DEEPCAD\_BIN**       executable directory
	+ **DEEPCAD\_DATA**    sample data folder



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. 
  
Sample session on a GPU node:



```

[user@biowulf ~]$ **sinteractive --mem=48g -c4 --gres=gpu:k80:1,lscratch:10**
[user@cn4219 ~]$ **module load deepcad**
[+] Loading singularity  3.8.4  on cn4219
[+] Loading cuDNN 7.0  libraries...
[+] Loading CUDA Toolkit  9.0.176 ...
[+] Loading deepcad  20210826

```

Basic usage:

```

[user@cn4219 ~]$  **main.py -h**
...
usage: main.py [-h] [--img_h IMG_H] [--img_w IMG_W] [--img_s IMG_S]
               [--img_c IMG_C] [--gap_h GAP_H] [--gap_w GAP_W] [--gap_s GAP_S]
               [--normalize_factor NORMALIZE_FACTOR]
               [--train_epochs TRAIN_EPOCHS] [--datasets_path DATASETS_PATH]
               [--datasets_folder DATASETS_FOLDER]
               [--DeepCAD_model_folder DEEPCAD_MODEL_FOLDER]
               [--results_folder RESULTS_FOLDER] [--GPU GPU]
               [--is_training IS_TRAINING] [--lr LR]
               [--train_datasets_size TRAIN_DATASETS_SIZE]
               [--select_img_num SELECT_IMG_NUM]

optional arguments:
  -h, --help            show this help message and exit
  --img_h IMG_H         the height of patch stack
  --img_w IMG_W         the width of patch stack
  --img_s IMG_S         the image number of patch stack
  --img_c IMG_C         the channel of image
  --gap_h GAP_H         the height of patch gap
  --gap_w GAP_W         the width of patch gap
  --gap_s GAP_S         the image number of patch gap
  --normalize_factor NORMALIZE_FACTOR
                        Image normalization factor
  --train_epochs TRAIN_EPOCHS
                        train epochs
  --datasets_path DATASETS_PATH
                        the name of your project
  --datasets_folder DATASETS_FOLDER
                        the folders for datasets
  --DeepCAD_model_folder DEEPCAD_MODEL_FOLDER
                        the folders for DeepCAD(pb) model
  --results_folder RESULTS_FOLDER
                        the folders for results
  --GPU GPU             the index of GPU you will use for computation
  --is_training IS_TRAINING
                        train or test
  --lr LR               initial learning rate
  --train_datasets_size TRAIN_DATASETS_SIZE
                        actions: train or predict
  --select_img_num SELECT_IMG_NUM
                        actions: train or predict

```

Create an input data folder and run DeepCAD on the data in that folder:

```

[user@cn4219 ~]$  **mkdir -p datasets/my\_data**
[user@cn4219 ~]$  **cp $DEEPCAD\_DATA/\* datasets/my\_data**
[user@cn4219 ~]$  **main.py --GPU 0 --img\_h 64 --img\_w 64 --img\_s 320 --train\_epochs 30 --datasets\_folder my\_data --normalize\_factor 1 --lr 0.00005 --train\_datasets\_size 100 &**
list(os.walk(im_folder, topdown=False)) ----->  [('datasets//my_data', [], ['ForTestPlugin_256x256x1000.tif'])]
stack_num ----->  1

TiffPage 0: TypeError: read_bytes() missing 3 required positional arguments: 'dtype', 'count', and 'offsetsize'
whole_w ----->  256
whole_h ----->  256
whole_s ----->  1000
w_num ----->  4
h_num ----->  4
s_num ----->  7
gap_s ----->  60
2021-11-19 15:37:52.507610: I tensorflow/core/platform/cpu_feature_guard.cc:137] Your CPU supports instructions that this TensorFlow binary was not compiled to use: SSE4.1 SSE4.2 AVX AVX2 FMA
2021-11-19 15:37:52.507957: E tensorflow/stream_executor/cuda/cuda_driver.cc:406] failed call to cuInit: CUresult(-1)
2021-11-19 15:37:52.508015: I tensorflow/stream_executor/cuda/cuda_diagnostics.cc:158] retrieving CUDA diagnostic information for host: cn4219
2021-11-19 15:37:52.508042: I tensorflow/stream_executor/cuda/cuda_diagnostics.cc:165] hostname: cn4219
2021-11-19 15:37:52.508197: I tensorflow/stream_executor/cuda/cuda_diagnostics.cc:189] libcuda reported version is: Not found: was unable to find libcuda.so DSO loaded into this program
2021-11-19 15:37:52.508245: I tensorflow/stream_executor/cuda/cuda_diagnostics.cc:369] driver version file contents: """NVRM version: NVIDIA UNIX x86_64 Kernel Module  470.82.01  Wed Oct 27 21:21:55 UTC 2021
GCC version:  gcc version 4.8.5 20150623 (Red Hat 4.8.5-28) (GCC)
"""
2021-11-19 15:37:52.508289: I tensorflow/stream_executor/cuda/cuda_diagnostics.cc:193] kernel reported version is: 470.82.1
--- Epoch  0  --- Step  0 / 112  --- L1_loss  382.85928  --- L2_loss  362411.6  --- Time  1013.0063416957855
train_input --->  14906.737 --->  0.0
output_img --->  1.1468295 --->  -0.42850912
/usr/local/apps/deepcad/20210826/src/main.py:89: UserWarning: results//unet3d_my_data_20211119-1537//0_0_my_data_ForTestPlugin_256x256x1000_x1_y2_z3_output.tif is a low contrast image
  io.imsave(result_name, output_img.transpose(2,0,1))
/usr/local/apps/deepcad/20210826/src/main.py:90: UserWarning: results//unet3d_my_data_20211119-1537//0_0_my_data_ForTestPlugin_256x256x1000_x1_y2_z3_noise1.tif is a low contrast image
  io.imsave(noise_img1_name, train_input.transpose(2,0,1))
/usr/local/apps/deepcad/20210826/src/main.py:91: UserWarning: results//unet3d_my_data_20211119-1537//0_0_my_data_ForTestPlugin_256x256x1000_x1_y2_z3_noise2.tif is a low contrast image
  io.imsave(noise_img2_name, train_GT.transpose(2,0,1))
...
[user@cn4219 ~]$  **nvidia-smi**
Fri Nov 19 15:37:53 2021
+-----------------------------------------------------------------------------+
| NVIDIA-SMI 470.82.01    Driver Version: 470.82.01    CUDA Version: 11.4     |
|-------------------------------+----------------------+----------------------+
| GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
| Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
|                               |                      |               MIG M. |
|===============================+======================+======================|
|   0  Tesla K80           On   | 00000000:8A:00.0 Off |                    0 |
| N/A   79C    P0   140W / 149W |  10955MiB / 11441MiB |    100%      Default |
|                               |                      |                  N/A |
+-------------------------------+----------------------+----------------------+

+-----------------------------------------------------------------------------+
| Processes:                                                                  |
|  GPU   GI   CI        PID   Type   Process name                  GPU Memory |
|        ID   ID                                                   Usage      |
|=============================================================================|
|    0   N/A  N/A     21956      C   ...a/envs/deepcad/bin/python    10950MiB |
+-----------------------------------------------------------------------------+
...
--- Epoch  1  --- Step  0 / 112  --- L1_loss  335.09457  --- L2_loss  224135.25  --- Time  1144.084424495697
train_input --->  6691.333 --->  0.0
output_img --->  1.9556417 --->  0.14389074
/usr/local/apps/deepcad/20210826/src/main.py:89: UserWarning: results//unet3d_my_data_20211120-0843//1_0_my_data_ForTestPlugin_256x256x1000_x0_y3_z4_output.tif is a low contrast image
  io.imsave(result_name, output_img.transpose(2,0,1))
/usr/local/apps/deepcad/20210826/src/main.py:90: UserWarning: results//unet3d_my_data_20211120-0843//1_0_my_data_ForTestPlugin_256x256x1000_x0_y3_z4_noise1.tif is a low contrast image
  io.imsave(noise_img1_name, train_input.transpose(2,0,1))
/usr/local/apps/deepcad/20210826/src/main.py:91: UserWarning: results//unet3d_my_data_20211120-0843//1_0_my_data_ForTestPlugin_256x256x1000_x0_y3_z4_noise2.tif is a low contrast image
  io.imsave(noise_img2_name, train_GT.transpose(2,0,1))
...
--- Epoch  20  --- Step  0 / 112  --- L1_loss  281.3424  --- L2_loss  174512.05  --- Time  21978.92845606804
train_input --->  7301.258 --->  0.0
output_img --->  25.690779 --->  6.511296
/usr/local/apps/deepcad/20210826/src/main.py:89: UserWarning: results//unet3d_my_data_20211120-0843//20_0_my_data_ForTestPlugin_256x256x1000_x1_y0_z3_output.tif is a low contrast image
  io.imsave(result_name, output_img.transpose(2,0,1))
/usr/local/apps/deepcad/20210826/src/main.py:90: UserWarning: results//unet3d_my_data_20211120-0843//20_0_my_data_ForTestPlugin_256x256x1000_x1_y0_z3_noise1.tif is a low contrast image
  io.imsave(noise_img1_name, train_input.transpose(2,0,1))
/usr/local/apps/deepcad/20210826/src/main.py:91: UserWarning: results//unet3d_my_data_20211120-0843//20_0_my_data_ForTestPlugin_256x256x1000_x1_y0_z3_noise2.tif is a low contrast image
  io.imsave(noise_img2_name, train_GT.transpose(2,0,1))
...

```

Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. deepcad.sh). For example:



```

#!/bin/bash
set -e
module load deepcad 
mddir datasets/my_data 
cp $DEEPCAD_DATA/* datasets/my_data
main.py --GPU 0 --img_h 64 --img_w 64 --img_s 320 --train_epochs 30 --datasets_folder my_data --normalize_factor 1 --lr 0.00005 --train_datasets_size 100

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
**sbatch deepcad.sh**
```





