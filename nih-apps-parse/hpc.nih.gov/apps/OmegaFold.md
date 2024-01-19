

document.querySelector('title').textContent = 'OmegaFold: High-resolution de novo Structure Prediction from Primary Sequence ';
**OmegaFold: High-resolution de novo Structure Prediction from Primary Sequence** 


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



OmegaFold is the first computational method to successfully predict high-resolution protein structure 
from a single primary sequence alone. 
Using a new combination of a protein language model 
that allows us to make predictions from single sequences 
and a geometry-inspired transformer model trained on protein structures, 
OmegaFold outperforms RoseTTAFold 
and achieves similar prediction accuracy to AlphaFold2 on recently released structures. 



### References:


* Ruidong Wua, Fan Dinga, Rui Wanga, Rui Shena, Xiwen Zhanga, Shitong Luoa, Chenpeng Sua, Zuofan Wua, Qi Xieb,
Bonnie Bergerc, Jianzhu Maa and Jian Penga   

*High-resolution de novo structure prediction from primary sequence.*   

[bioRxiv preprint doi: https://doi.org/10.1101/2022.07.21.500999; posted July 22, 2022. T](https://www.biorxiv.org/content/10.1101/2022.07.21.500999v1.abstract)


Documentation
* [OmegaFold GitHub page](https://github.com/HeliXonProtein/OmegaFold)


Important Notes
* Module Name: OmegasFold (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **OMEGAFOLD\_HOME**  installation directory
	+ **OMEGAFOLD\_BIN**   executable directory
	+ **OMEGAFOLD\_DATA**  test data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
  

The OmegaFold application installed on Biowulf is supposed to be run on a GPU node. 
As the first step, please allocate an [interactive session](/docs/userguide.html#int):



```

[user@biowulf]$ **sinteractive --gres=gpu:a100:1,lscratch:10 --mem=100g -c16**
[user@cn1113 ~]$ **module load OmegaFold** 
[+] Loading singularity  3.8.5-1  on cn1113
[+] Loading CUDA Toolkit  11.3.0  ...
[+] Loading cuDNN/8.2.1/CUDA-11.3 libraries...
[+] Loading OmegaFold  1.1.0
[user@cn1113 user]$ **omegafold --help**
usage: omegafold [-h] [--num_cycle NUM_CYCLE] [--subbatch_size SUBBATCH_SIZE] [--device DEVICE]
                 [--weights_file WEIGHTS_FILE] [--weights WEIGHTS]
                 [--pseudo_msa_mask_rate PSEUDO_MSA_MASK_RATE] [--num_pseudo_msa NUM_PSEUDO_MSA]
                 [--allow_tf32 ALLOW_TF32]
                 input_file output_dir

Launch OmegaFold and perform inference on the data. Some examples (both the input and output
files) are included in the Examples folder, where each folder contains the output of each
available model from model1 to model3. All of the results are obtained by issuing the general
command with only model number chosen (1-3).

positional arguments:
  input_file            The input fasta file
  output_dir            The output directory to write the output pdb files. If the directory does
                        not exist, we just create it. The output file name follows its unique
                        identifier in the rows of the input fasta file"

optional arguments:
  -h, --help            show this help message and exit
  --num_cycle NUM_CYCLE
                        The number of cycles for optimization, default to 10
  --subbatch_size SUBBATCH_SIZE
                        The subbatching number, the smaller, the slower, the less GRAM
                        requirements. Default is the entire length of the sequence. This one takes
                        priority over the automatically determined one for the sequences
  --device DEVICE       The device on which the model will be running, default to the accelerator
                        that we can find
  --weights_file WEIGHTS_FILE
                        The model cache to run
  --weights WEIGHTS     The url to the weights of the model
  --pseudo_msa_mask_rate PSEUDO_MSA_MASK_RATE
                        The masking rate for generating pseudo MSAs
  --num_pseudo_msa NUM_PSEUDO_MSA
                        The number of pseudo MSAs
  --allow_tf32 ALLOW_TF32
                        if allow tf32 for speed if available, default to True

```

Copy testing data to your current folder: 

```

[user@cn1113 user]$ **cp $OMEGAFOLD\_DATA/\* .**

```

Run the omegafold executable on testingf data:

```

[user@cn1113 user]$ **omegafold pi3k.fa outdir &** 
[1] 530515
[user@cn1113 OmegaFold]$ INFO:root:Loading weights from /home/user/.cache/omegafold_ckpt/model.pt
INFO:root:Constructing OmegaFold
INFO:root:Reading pi3k.fa
INFO:root:Predicting 1th chain in pi3k.fa
INFO:root:724 residues in this chain.
...
[user@cn1113 OmegaFold]$ nvidia-smi
+-----------------------------------------------------------------------------+
| NVIDIA-SMI 470.82.01    Driver Version: 470.82.01    CUDA Version: 11.4     |
|-------------------------------+----------------------+----------------------+
| GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
| Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
|                               |                      |               MIG M. |
|===============================+======================+======================|
|   0  NVIDIA A100-SXM...  On   | 00000000:C7:00.0 Off |                    0 |
| N/A   61C    P0   367W / 400W |  23136MiB / 81251MiB |    100%      Default |
|                               |                      |             Disabled |
+-------------------------------+----------------------+----------------------+

+-----------------------------------------------------------------------------+
| Processes:                                                                  |
|  GPU   GI   CI        PID   Type   Process name                  GPU Memory |
|        ID   ID                                                   Usage      |
|=============================================================================|
|    0   N/A  N/A    530539      C   /usr/bin/python3                23133MiB |
+-----------------------------------------------------------------------------+

...
INFO:root:Finished prediction in 185.76 seconds.
INFO:root:Saving prediction to outdir/sp|P27986|P85A_HUMAN Phosphatidylinositol 3-kinase regulatory subunit alpha OS=Homo sapiens OX=9606 GN=PIK3R1 PE=1 SV=2.pdb
INFO:root:Saved
INFO:root:Predicting 2th chain in pi3k.fa
INFO:root:1068 residues in this chain.
INFO:root:Finished prediction in 529.51 seconds.
INFO:root:Saving prediction to outdir/sp|P42336|PK3CA_HUMAN Phosphatidylinositol 4,5-bisphosphate 3-kinase catalytic subunit alpha isoform OS=Homo sapiens OX=9606 GN=PIK3CA PE=1 SV=2.pdb
INFO:root:Saved
INFO:root:Done!

[user@cn1113 OmegaFold]$ **ls outdir**
'sp|P27986|P85A_HUMAN Phosphatidylinositol 3-kinase regulatory subunit alpha OS=Homo sapiens OX=9606 GN=PIK3R1 PE=1 SV=2.pdb'
'sp|P42336|PK3CA_HUMAN Phosphatidylinositol 4,5-bisphosphate 3-kinase catalytic subunit alpha isoform OS=Homo sapiens OX=9606 GN=PIK3CA PE=1 SV=2.pdb'

```

En the interactive session:

```

[user@cn1113 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226

```





