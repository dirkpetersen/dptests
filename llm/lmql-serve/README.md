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

Before we start we need an interactive session on an HPC node, let's request 4 CPU cores, 64GB memory and a single GPU for one day

```
srun -c 4 --mem 64G -p gpu --gpus 1 -t 1-0 --pty bash
```

Then we need to install a conda environment or pixi project with the right software dependencies. We need to make sure that `lmql` can use GPUs locally, for which the llama-cpp-python package must be compiled to support GPUs. lmql requires at least Python verison 3.10 while llama-cpp-python supports GCC version 13 or older.

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

Now we create a slurm batch job that loads 2 inference servers, one with a specific model and the other one in auto mode so that you can load a model. 

vi lmql-server.sub 

```
#! /bin/bash
#SBATCH --job-name "lmqlsvr"
#SBATCH --time 1-00:00:00
#SBATCH --cpus-per-task 4
#SBATCH --mem 64G
#SBATCH --partition gpu 
#SBATCH --gpus a40:1 
#SBATCH --output lmqlsvr-%J.out
#SBATCH --error lmqlsvr-%J.err

nvidia-smi
host=$(hostname)
port=$(/usr/local/bin/get-open-port.py)

pixi shell --manifest-path /my/conda/environments/lmqlsvr/pixi.toml

lmql serve-model --cuda --host ${host} --port ${port} llama.cpp:/arc/scratch1/dpcri/llm/mistral-7b-instruct-v0.1.Q4_K_M.gguf

# or you can use auto mode if you want to load a model at runtime 
# lmql serve-model --cuda --host ${host} --port ${port} auto

```

then you run the job and check the job output file to see which node and port the server is listening on:

```
sbatch lmql-server.sub
tail lmqlsvr-1234567.out

[Serving LMTP endpoint on ws://cnode-11-16:19000/]
[Loading llama.cpp model from llama.cpp:/arc/scratch1/dpcri/llm/mistral-7b-instruct-v0.1.Q4_K_M.gguf  with  {'device_map': 'auto'} ]
```

Now you can access this model from the login node, for example by running an interactive session using VS code connecting your lmql client to `cnode-11-16:19000` or whatever the node port combination is.

**Warning: this is a security hole as anyone on the cluster will be able to access the port if they know the port number and node and they will have access to your models if they know lmql.**


