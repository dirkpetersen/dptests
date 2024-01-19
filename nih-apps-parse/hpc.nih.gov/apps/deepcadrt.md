

document.querySelector('title').textContent = 'DeepCAD: Deep self-supervised learning for calcium imaging denoising..';
**DeepCAD-RT: Real-time denoising of fluorescence time-lapse imaging using deep self-supervised learning.**


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



DeepCAD-RT can denoise fluorescence time-lapse images with rapid processing speed that can be incorporated with the microscope acquisition system to achieve real-time denoising. Our method is based on deep self-supervised learning and the original low-SNR data can be directly used for training convolutional networks, making it particularly advantageous in functional imaging where the sample is undergoing fast dynamics and capturing ground-truth data is hard or impossible. We have demonstrated extensive experiments including calcium imaging in mice, zebrafish, and flies, cell migration observations, and the imaging of a new genetically encoded ATP sensor, covering both 2D single-plane imaging and 3D volumetric imaging. 



### References:


 Li X, Li Y, Zhou Y, Wu J, Zhao Z, Fan J, Deng F, Wu Z, Xiao G, He J, Zhang Y, Zhang G, Hu X, Chen X, Zhang Y, Qiao H, Xie H, Li Y, Wang H, Fang L, Dai Q. *Real-time denoising enables high-sensitivity fluorescence time-lapse imaging beyond the shot-noise limit.* Nat Biotechnol. 2023 Feb;41(2):282-292. doi: 10.1038/s41587-022-01450-8. Epub 2022 Sep 26. PMID: 36163547; PMCID: PMC9931589.[Journal](https://pubmed.ncbi.nlm.nih.gov/36163547/) 

Documentation
* [DeepCAD-RT Github page](https://github.com/cabooster/DeepCAD-RT)


Common pitfalls
[top](#top)


**cv2.error**
There is no GUI installed, so if you see some errors regardinding "cv2.error: OpenCV(4.7.0) /io/opencv/modules/highgui/src/window.cpp:1255: error: (-2:Unspecified error) The function is not implemented. Rebuild the library with Windows, GTK+ 2.x or Cocoa support. If you are on Ubuntu or Debian, install libgtk2.0-dev and pkg-config, then re-run cmake or configure script in function 'cvNamedWindow'", please modify the code to set visulization off.
 
```

		     visualize_images_per_epoch = False
		     display_images = False
	     
```


Important Notes
* Module Name: deepcadrt (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* This is just the command line venison with Jupyter, the GUI interface and Matlab GUI is not included.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. 
  
Sample session on a GPU node:



```

[user@biowulf ~]$ **sinteractive --mem=48g -c4 --gres=gpu:v100:1,lscratch:10**
[user@cn4219 ~]$ **module load deepcadrt**
Loading deepcadrt  0.1.0

```

Basic usage:

```

[user@cn4219 ~]$  **python-deepcadrt**

from deepcad.train_collection import training_class
from deepcad.movie_display import display
from deepcad.utils import get_first_filename,download_demo


```

Download the data from DeepCAD-RT git repo:

```

[user@cn4219 ~]$  **mkdir -p /data/$USER/deepcadrt\_data; cd /data/$USER/deepcadrt\_data**
[user@cn4219 ~]$  **git clone https://github.com/cabooster/DeepCAD-RT**
[user@cn4219 ~]$  **cd DeepCAD-RT/DeepCAD\_RT\_pytorch/**
[user@cn4219 ~]$  **python-deepcadrt demo\_train\_pipeline.py**
python-deepcadrt demo_train_pipeline.py
Training parameters ----->
{'overlap_factor': 0.25, 'datasets_path': 'datasets/fish_localbrain_demo', 'n_epochs': 10, 'fmap': 16, 'output_dir': './results', 'pth_dir': './pth', 'onnx_dir': './onnx', 'batch_size': 1, 'patch_t': 150, 'patch_x': 150, 'patch_y': 150, 'gap_y': 112, 'gap_x': 112, 'gap_t': 112, 'lr': 5e-05, 'b1': 0.5, 'b2': 0.999, 'GPU': '0', 'ngpu': 1, 'num_workers': 4, 'scale_factor': 1, 'train_datasets_size': 6000, 'select_img_num': 100000, 'test_datasize': 400, 'visualize_images_per_epoch': True, 'save_test_images_per_epoch': True, 'colab_display': False, 'result_display': ''}
Image list for training ----->
Total stack number ----->  1
Noise image name ----->  fish_localbrain.tif
...

```

Run with jupyter notebook:
lease set up a tunnel as the jupyter webpage: [Jupyter](https://hpc.nih.gov/apps/jupyter.html)
Then run the same command without loading the jupyter module:

```

jupyter notebook --ip localhost --port $PORT1 --no-browser

```

Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. deepcadrt.sh). For example:



```

#!/bin/bash
set -e
module load deepcadrt  
python-deepcadrt demo_train_pipeline.py

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
**sbatch --gres=gpu:v100:1,lscratch:20 deepcadrt.sh** 
```





