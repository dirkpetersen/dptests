

document.querySelector('title').textContent = 'DeepMM: full-length protein structure determination from cryo-EM maps using deep learning.';
DeepMM: full-length protein structure determination   
 from cryo-EM maps using deep learning.


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



DeepMM implements fully automated de novo structure modeling method, MAINMAST, 
which builds three-dimensional models of a protein from a near-atomic resolution EM map. 
The method directly traces the protein’s main-chain and identifies Cα
positions as tree-graph structures in the EM map.




DeepMM uses MAINMAST to trace the main chain paths on main chain probability maps. 

### References:


* Terashi, Genki, and Daisuke Kihara   

 *De novo main-chain modeling for EM maps using MAINMAST.*  

[Nature communications](https://www.nature.com/articles/s41467-018-04053-7)  (2018), **9**(1), 1618. 
* Terashi, Genki, and Daisuke Kihara  

 *De novo main-chain modeling with MAINMAST in 2015/2016 EM Model Challenge.*   

[Journal of Structural Biology](https://www.sciencedirect.com/science/article/pii/S1047847718301710)  (2018), **204**(2), 351-359.


Documentation
* [DeepMM Github page](https://github.com/JiahuaHe/DeepMM)


Important Notes
* Module Name: DeepMM (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **DEEPMM\_HOME**  installation directory
	+ **DEEPMM\_BIN**       executable directory
	+ **DEEPMM\_SRC**     a folder containing the source code
	+ **DEEPMM\_DATA**    a folder containing sample data



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Before running an interactive session, install Jackal software locally in your account:   

1) following the link:   
 https://honiglab.c2b2.columbia.edu/software/cgi-bin/software.pl?input=Jackal   

 manually download the file jackal\_64bit.tar.gz to your local computer   

2) sftp/transfer the file to Biowulf   

3) ungzip/untar it:   

 tar -zxf jackal\_64bit.tar.gz   
 
 and   

4) move the resulting folder jackal\_64bit   
 to the directory: /data/$USER 


Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=8g --gres=gpu:k80:1,lscratch:10 -c4** 
[user@cn4213 ~]$ **module load DeepMM** 
[+] Loading CUDA Toolkit  10.2.89  ...
[+] Loading cuDNN/7.6.5/CUDA-10.2 libraries...
[+] Loading blast 2.11.0+  ...
[+] Loading DeepMM  20210722

```

Set the environment variables:

```

[user@cn4213 ~]$  **export JACKALDIR=/data/$USER/jackal\_64bit**
[user@cn4213 ~]$  **PATH=/data/$USER/jackal\_64bit/bin:$PATH** 
[user@cn4213 ~]$  **git clone https://github.com/JiahuaHe/DeepMM** 

```

Run DeepMM on the sample data:

```

[user@cn4213 ~]$ **cd DeepMM/5185**
[user@cn4213 ~]$ **preprocess.py -m 5185\_zoned.mrc -n map.npz -t 3.21**  
[INFO] 10733 voxels with density value greater than 3.21000 are retained.
[INFO] 10733 voxels are saved to map.npz

[user@cn4213 ~]$ **pred\_MCCA.py -n map.npz -m predMC.mrc -c predCA.mrc** 
LOAD] loading data...
[PREDICT] predicting...
100%|████████████████████████████████████████████████████████████████████████████████████████████████| 84/84 [00:04<00:00, 17.43it/s]

[user@cn4213 ~]$  **mrc2situs.py -m 5185\_zoned.mrc -s map.situs** 
[user@cn4213 ~]$  **mrc2situs.py -m predMC.mrc -s predMC.situs**   
[user@cn4213 ~]$  **mrc2situs.py -m predCA.mrc -s predCA.situs**   

[user@cn4213 ~]$  **getldp predMC.situs predCA.situs > LDP.pdb**   

[user@cn4213 ~]$  **getvox map.situs LDP.pdb > LDP.mcv** 

[user@cn4213 ~]$  **pred\_AASS.py -m LDP.mcv -l LDP.pdb > LDP\_AASS.pdb** 
100%|█████████████████████████████████████████████████████████████████████| 7/7 [00:03<00:00,  2.18it/s]

[user@cn4213 ~]$  **trace LDP\_AASS.pdb -nrd 100 > paths.pdb** 

[user@cn4213 ~]$  **align paths.pdb seq.fasta.spd3** 
# ALIGN.F
 # Align sequence to main-chain paths
 # Parameter setting:
 # wca =   1.600000
 # whelix =   1.000000
 # wsheet =  0.7000000
 # wcoil =  0.8000000
 # waa =   1.000000
 # wss =  0.5000000
 # rmsd =   5.000000
 # nout =          10
 # prefix = model
 # Reading LDP paths from file paths.pdb
 # Path            1  consists of          565  LDPs.
 # Path            2  consists of          560  LDPs.
 # Path            3  consists of          534  LDPs.
 # Path            4  consists of          538  LDPs.
 # Path            5  consists of          511  LDPs.
 # Path            6  consists of          520  LDPs.
 # Path            7  consists of          543  LDPs.
 # Path            8  consists of          484  LDPs.
 # Path            9  consists of          551  LDPs.
 # Path           10  consists of          475  LDPs.
 # Number of paths:           10
 # Reading SS prediction from file seq.fasta.spd3
 # Reading predicted SS from .SPD3 file seq.fasta.spd3
 # Number of residues:          155
 # Start alignment!
 # Model    1   62.156
 # Model    2   57.567
 # Model    3   78.881
 # Model    4   76.093
 # Model    5   78.612
 # Model    6   90.198
 # Model    7   87.193
 # Model    8   73.144
 # Model    9   71.957
 # Model   10  100.461
 # Model   11  112.469
 # Model   12  108.363
 # Model   13  139.715
 # Model   14  130.736
 # Model   15  123.928
 # Model   16  116.492
 # Model   17   34.216
 # Model   18   35.854
 # Model   19   27.489
 # Model   20   44.658
 # Model   21   45.993
 # Model   22   52.989
 # Model   23   42.962
 # Model   24   34.723
 # Model   25  110.321
... 
 # Model  158   59.206
 # Model  159   40.211
 # Model  160   27.201
 # write models into one .pdb file

```

End the interactive session:

```

[user@cn4213 ~]$ **exit**
[user@biowulf ~]$

```





