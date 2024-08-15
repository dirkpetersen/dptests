# Running an lmql inference server on slurm GPU

The instructions work on OHSU ARC cluster 
Last change: 8/15/2024 

## prepare conda pixi / envionment 

Warning: Conda environments in AI/ML research can be massive, you should 1. move your application `~/.cache` folder away from your home directory to a shared scratch space and 2. also create your conda envionments in a fast and large file system. On arc do this: 

```
mkdir /arc/scratch1/yoursharedfolder/.cache
rclone copy ~/.cache /arc/scratch1/yoursharedfolder/.cache --progress --metadata --links --no-update-modtime
rm -rf ~/.cache
ln -s /arc/scratch1/yoursharedfolder/.cache ~/.cache
```

Before we start, we need an interactive session on an HPC node. Let's request 4 CPU cores, 64GB memory and a single GPU for one day

```
srun -c 4 --mem 64G -p gpu --gpus 1 -t 1-0 --pty bash
```

Then we need to install a conda environment or pixi project with the right software dependencies. We need to make sure that `lmql` can use GPUs locally, for which the llama-cpp-python package must be compiled to support GPUs. lmql requires at least Python verison 3.10 while llama-cpp-python supports up to GCC version 13.

```
curl -fsSL https://pixi.sh/install.sh | bash
```

Once you can use the pixi command (you might have to login again)  you run 

```
cd /my/conda/environments
pixi init lmqlsvr
cd lmqlsvr 
pixi add python=3.12 pip ipython gcc=13 gxx=13 cuda=12.5
pixi add --pypi lmql[hf] transformers
pixi shell
```

then you can either install llama-cpp-python from scratch :

```
CMAKE_ARGS="-DGGML_CUDA=on" pip install llama-cpp-python
```

or you install a precompiled wheel on ARC

```
pip install /arc/scratch1/dpcri/llm/whl/cuda-12.5/llama_cpp_python-0.2.88-cp312-cp312-linux_x86_64.whl
```

Now we create a slurm batch job that loads an inference server, you can either load one with a specific model or you choose auto mode so that you can load different models at run time. 

vi lmql-server.sub 

```
#! /bin/bash
#SBATCH --job-name "lmqlsvr"
#SBATCH --time 1-00:00:00
#SBATCH --cpus-per-task 4
#SBATCH --mem 64G
#SBATCH --partition gpu 
#SBATCH --gpus a40:1 
#SBATCH --output lmqlsvr.out
#SBATCH --error lmqlsvr.err

# nvidia-smi # for debugging
host=$(hostname)
port=$(/usr/local/bin/get-open-port.py)

pixi run --manifest-path /my/conda/environments/lmqlsvr/pixi.toml \
 lmql serve-model --cuda --host ${host} --port ${port} \
 llama.cpp:/arc/scratch1/dpcri/llm/mistral-7b-instruct-v0.1.Q4_K_M.gguf

# or you can use auto mode if you want to load a model at runtime 
#pixi run --manifest-path /my/conda/environments/lmqlsvr/pixi.toml \
# lmql serve-model --cuda --host ${host} --port ${port} auto

```

then you run the job and check the job output file to see which node and port the server is listening on:

```
sbatch lmql-server.sub
tail lmqlsvr.out

[Serving LMTP endpoint on ws://cnode-11-16:19000/]
[Loading llama.cpp model from llama.cpp:/arc/scratch1/dpcri/llm/mistral-7b-instruct-v0.1.Q4_K_M.gguf  with  {'device_map': 'auto'} ]
```

Now you can access this model from the login node, for example by running an interactive session using VS code connecting your lmql client to `cnode-11-16:19000` or whatever the node port combination is.

**Warning: Not binding a service to localhost  is a security hole as anyone on the cluster will be able to access the port if they know the port number and node and they will have access to your models if they know lmql.**


To verify that it is in fact using a GPU you can directly login to the node with ssh for monitoring purposes since you have an active job running on that node. After logging in run `nvidia-smi`. You should see that one of your GPUs is using about 450 MB GPU memory with this specific model


```
ssh cnode-11-16
nvidia-smi

Thu Aug 15 15:50:07 2024
+-----------------------------------------------------------------------------------------+
| NVIDIA-SMI 555.42.06              Driver Version: 555.42.06      CUDA Version: 12.5     |
|-----------------------------------------+------------------------+----------------------+
| GPU  Name                 Persistence-M | Bus-Id          Disp.A | Volatile Uncorr. ECC |
| Fan  Temp   Perf          Pwr:Usage/Cap |           Memory-Usage | GPU-Util  Compute M. |
|                                         |                        |               MIG M. |
|=========================================+========================+======================|
|   2  NVIDIA A40                     On  |   00000000:41:00.0 Off |                    0 |
|  0%   37C    P0             74W /  300W |     451MiB /  46068MiB |      0%      Default |
|                                         |                        |                  N/A |
+-----------------------------------------+------------------------+----------------------+

+-----------------------------------------------------------------------------------------+
| Processes:                                                                              |
|  GPU   GI   CI        PID   Type   Process name                              GPU Memory |
|        ID   ID                                                               Usage      |
|=========================================================================================|
|    2   N/A  N/A   1743472      C   ...r/.pixi/envs/default/bin/python3.12          0MiB |
+-----------------------------------------------------------------------------------------+ 
```

in addition to output file the hostname and port information has also been saved in the Comment field of the Slurm job and you can retrieve it from there, assuming your slurm job is 123456.

```
scontrol show job 123456 | sed -n 's/.*Comment=\([^[:space:]]*\).*/\1/p'

lmqlsvr:cnode-11-16:19000
```


