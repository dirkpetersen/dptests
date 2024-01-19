

document.querySelector('title').textContent = 'B-SOiD: going from position to behaviors ';
**B-SOiD: going from position to behaviors** 


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



B-SOiD (Behavioral Segmentation in Deeplabcut) is an unsupervised learning
algorithm that serves to discover and classify behaviors that are not pre-defined by users.
It segregates statistically different, sub-second rodent behaviors with a single
bottom-up perspective video cameraR by performing a novel expectation maximization
fitting of Gaussian mixture models on t-Distributed Stochastic Neighbor Embedding (t-SNE).



### References:


* Alexander I. Hsu and Eric A. Yttri   

*B-SOiD, an open-source unsupervised algorithm for identification and fast prediction of behaviors.*   

[Nature Communications volume 12, Article number: 5188 (2021)](https://www.nature.com/articles/s41467-021-25420-x)


Documentation
* [B-SOiD GitHub page](https://github.com/YttriLab/B-SOID)
* [B-SOiD python tutorial](https://github.com/YttriLab/B-SOID/blob/master/docs/bsoid_umap_tutorial.md)


Important Notes
* Module Name: B-SOiD (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **BSOID\_HOME**  installation directory
	+ **BSOID\_BIN**    executables directory
	+ **BSOID\_SRC**    source directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --gres=gpu:p100:1,lscratch:10 --mem=20g -c14**
[user@cn4469 ~]$ **module load b-soid** 
[+] Loading singularity  3.8.5-1  on cn0897
[+] Loading b-soid  1.5.1

```

Create and enter the project folder, e.g.: 

```

[user@cn4469 user]$ **mkdir /data/$USER/B-SOID && cd /data/$USER/B-SOID** 

```

Create and populate data subfolders, e.g.

```

[user@cn4469 user]$ **mkdir -p datasets/Train1 datasets/Train2 datasets/Data1/bsoid\_umap\_beta** 
[user@cn4469 user]$ **cp $BSOID\_DATA/bl6c57/041919/\* datasets/Train1** 
[user@cn4469 user]$ **cp $BSOID\_DATA/bl6c57/042219/\* datasets/Train2** 

```

Copy the B-SOID/1.3 source code to your current directory

```

[user@cn4469 user]$ **cp -r $BSOID\_SRC/\* .**

```

Customoze the configuration files bsoid\_py/config/LOCAL\_CONFIG.py and bsoid\_umap/config/LOCAL\_CONFIG.py, e.g. by editing properly the lines containing the string "user".   
  

Now build your behaviaoral model:

```

[user@cn4469 ~]$ **ipython** 
In [1]: **from bsoid\_umap.config import \***

In [2]: **import bsoid\_umap.main**

In [3]: **f\_10fps, f\_10fps\_sc, umap\_embeddings, hdb\_assignments, soft\_assignments, soft\_clusters, nn\_classifier, scores, nn\_assignments = bsoid\_umap.main.build(TRAIN\_FOLDERS)**
2022-06-23 09:18:44 INFO     Importing CSV file 1 from folder 1
2022-06-23 09:18:45 INFO     Extracting likelihood value...
2022-06-23 09:18:45 INFO     Computing data threshold to forward fill any sub-threshold (x,y)...
100%|█████████████████████████████████████████████████████████████████| 6/6 [00:07<00:00,  1.32s/it]
2022-06-23 09:18:53 INFO     Done preprocessing (x,y) from file 1, folder 1.
2022-06-23 09:18:53 INFO     Importing CSV file 2 from folder 1
2022-06-23 09:18:54 INFO     Extracting likelihood value...
2022-06-23 09:18:54 INFO     Computing data threshold to forward fill any sub-threshold (x,y)...
100%|█████████████████████████████████████████████████████████████████| 6/6 [00:07<00:00,  1.33s/it]
2022-06-23 09:19:02 INFO     Done preprocessing (x,y) from file 2, folder 1.
2022-06-23 09:19:02 INFO     Processed 2 CSV files from folder: /Train1
2022-06-23 09:19:02 INFO     Importing CSV file 1 from folder 2
2022-06-23 09:19:04 INFO     Extracting likelihood value...
2022-06-23 09:19:04 INFO     Computing data threshold to forward fill any sub-threshold (x,y)...
100%|█████████████████████████████████████████████████████████████████| 6/6 [00:07<00:00,  1.32s/it]
2022-06-23 09:19:12 INFO     Done preprocessing (x,y) from file 1, folder 2.
2022-06-23 09:19:12 INFO     Importing CSV file 2 from folder 2
2022-06-23 09:19:13 INFO     Extracting likelihood value...
2022-06-23 09:19:13 INFO     Computing data threshold to forward fill any sub-threshold (x,y)...
100%|█████████████████████████████████████████████████████████████████| 6/6 [00:07<00:00,  1.31s/it]
2022-06-23 09:19:20 INFO     Done preprocessing (x,y) from file 2, folder 2.
2022-06-23 09:19:20 INFO     Processed 2 CSV files from folder: /Train2
2022-06-23 09:19:20 INFO     Processed a total of 4 CSV files, and compiled into a (4, 107999, 12) data list.
2022-06-23 09:19:21 INFO     Extracting features from CSV file 1...
2022-06-23 09:21:57 INFO     Extracting features from CSV file 2...
2022-06-23 09:24:35 INFO     Extracting features from CSV file 3...
2022-06-23 09:27:12 INFO     Extracting features from CSV file 4...
2022-06-23 09:29:50 INFO     Done extracting features from a total of 4 training CSV files.
2022-06-23 09:29:54 INFO     Done integrating features into 100ms bins from CSV file 1.
2022-06-23 09:29:59 INFO     Done integrating features into 100ms bins from CSV file 2.
2022-06-23 09:30:04 INFO     Done integrating features into 100ms bins from CSV file 3.
2022-06-23 09:30:09 INFO     Done integrating features into 100ms bins from CSV file 4.
2022-06-23 09:30:09 INFO     Transforming all 71996 instances from 36 D into 3 D
2022-06-23 09:36:03 INFO     Done non-linear transformation with UMAP from 36 D into 3 D.
2022-06-23 09:36:03 INFO     Running HDBSCAN on 71996 instances in 3 D space...
2022-06-23 09:36:07 INFO     Adjusting minimum cluster size to maximize cluster number...
2022-06-23 09:36:55 INFO     Done predicting labels for 71996 instances in 3 D space...
2022-06-23 09:36:56 INFO     Training feedforward neural network on randomly partitioned 80.0% of training data...
2022-06-23 10:18:13 INFO     Done training feedforward neural network mapping (57596, 36) features to (57596,) assignments.
/opt/conda/envs/b-soid/lib/python3.7/site-packages/sklearn/utils/deprecation.py:87: FutureWarning: Function plot_confusion_matrix is deprecated; Function `plot_confusion_matrix` is deprecated in 1.0 and will be removed in 1.2. Use one of the class methods: ConfusionMatrixDisplay.from_predictions or ConfusionMatrixDisplay.from_estimator.
  warnings.warn(msg, category=FutureWarning)
connect localhost port 6000: Connection refused
Non-normalized confusion matrix
[[14230     1]
 [    8   161]]
/opt/conda/envs/b-soid/lib/python3.7/site-packages/sklearn/utils/deprecation.py:87: FutureWarning: Function plot_confusion_matrix is deprecated; Function `plot_confusion_matrix` is deprecated in 1.0 and will be removed in 1.2. Use one of the class methods: ConfusionMatrixDisplay.from_predictions or ConfusionMatrixDisplay.from_estimator.
  warnings.warn(msg, category=FutureWarning)
Normalized confusion matrix
[[1.00e+00 7.03e-05]
 [4.73e-02 9.53e-01]]
/gpfs/gsfs7/users/user/B-SOID/bsoid_umap/train.py:216: UserWarning: Matplotlib is currently using agg, which is a non-GUI backend, so cannot show the figure.
  plt.show()
2022-06-23 10:18:29 INFO     Scored cross-validated feedforward neural network performance.
/gpfs/gsfs7/users/user/B-SOID/bsoid_umap/utils/visuals.py:44: UserWarning: Matplotlib is currently using agg, which is a non-GUI backend, so cannot show the figure.
  plt.show()
/gpfs/gsfs7/users/user/B-SOID/bsoid_umap/utils/visuals.py:60: UserWarning: Matplotlib is currently using agg, which is a non-GUI backend, so cannot show the figure.
  plt.show()
2022-06-23 10:18:42 INFO     Saved.

In [4]: **data\_new, fs\_labels = bsoid\_umap.main.run(PREDICT\_FOLDERS)**

```

The output results will be stored automatically in the subfolder Output.

```

**quit**
[user@cn4469 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226

```





