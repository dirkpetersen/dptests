# Running a llama-cpp-python inference server on HPC GPU

The instructions work on OHSU ARC cluster using Slurm
Last change: 8/17/2024 

## prepare software environment 

Warning: Software packages in AI/ML research can be massive, you should 1. move your application `~/.cache` folder away from your home directory to a shared scratch space and 2. also create your conda or venv environments in a fast and large file system. On ARC do this: 

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

confirm that you have a GPU available

```
nvidia-smi
```

## Option1: compile llama-cpp-python from scratch

To compile we need to install a conda environment or pixi project with the right software dependencies. We need to make sure that the llama-cpp-python package is compiled to use GPUs.  llama-cpp-python supports up to GCC version 13.

if you don't have pixi installed, run :

```
curl -fsSL https://pixi.sh/install.sh | bash
```

Once you can use the pixi command (you might have to login again)  you run 

```bash
cd /my/conda/environments
pixi init llama-cpp
cd llama-cpp
pixi add python=3.9 pip ipython gcc=13 gxx=13 cuda=12.6
pixi add --pypi openai
pixi shell
```

then you can either install llama-cpp-python from scratch like this

```
CMAKE_ARGS="-DGGML_CUDA=on" pip install --upgrade llama-cpp-python[server]
```

## Option 2: use the precompiled wheel 

or you install a precompiled llama-cpp-python wheel on ARC using a Python 3.9 package either in a venv or in the default path in your home directory. Python 3.9 is the default Python on ARC.  

```
python3 -m pip install --upgrade /arc/scratch1/dpcri/llm/whl/cuda-12.6/llama_cpp_python-0.2.89-cp39-cp39-linux_x86_64.whl[server]
python3 -m pip install openai
```

## Create a Slurm job 

Now we create a slurm batch job that loads an openai API compatible inference server with a specific model 

vi llama-cpp-server.sub 

```bash 
#! /bin/bash
#SBATCH --job-name "llama-cpp-server"
#SBATCH --time 1-00:00:00
#SBATCH --cpus-per-task 4
#SBATCH --mem 64G
#SBATCH --partition gpu 
#SBATCH --gpus a40:1 
#SBATCH --output llama-cpp-server.out

OPENAI_API_KEY=${OPENAI_API_KEY:-$(whoami)} # if empty use username
# nvidia-smi # for debugging
host=$(hostname)
port=$(/usr/local/bin/get-open-port.py)
scontrol update JobId=${SLURM_JOB_ID} Comment="llama-cpp-server:${host}:${port}"

python3 -m llama_cpp.server --api_key ${OPENAI_API_KEY} --cache true \
 --host ${host} --port ${port} --n_gpu_layers -1 --verbose false \
 --model Meta-Llama-3.1-70B-Instruct.Q4_K_M.gguf


```

if you have or require  conda / pixi environment, the last command in the submission script is a bit different: 

```
pixi run --manifest-path /my/conda/environments/llama-cpp/pixi.toml \
 HOST=${host} python3 -m llama_cpp.server --api_key ${OPENAI_API_KEY} --cache true \
 --host ${host} --port ${port} --n_gpu_layers -1 --verbose false \
 --model Meta-Llama-3.1-70B-Instruct.Q4_K_M.gguf   # This is the largest model fitting into a A40 GPU
```

`--n_gpu_layers -1` means that everything should be offloaded to GPU, use `--n_gpu_layers 0` to run all inference in CPU. Use a positive integer number to run partially in GPU and CPU. The higher the number, the more GPU 

## Run and monitor slurm job 

then you run the job and check the job output file to see which node and port the server is listening on:

```
sbatch llama-cpp-server.sub
tail llama-cpp-server.out

INFO:     Started server process [2422952]
INFO:     Waiting for application startup.
INFO:     Application startup complete.
INFO:     Uvicorn running on http://cnode-12-2:33775 (Press CTRL+C to quit)
```

To verify that it is in fact using a GPU, you can directly login to the node with ssh for monitoring purposes since you have an active job running on that node. After logging in run `nvidia-smi`. You should see that one of your GPUs is using about 450 MB GPU memory with this specific model


```
ssh cnode-12-2
nvidia-smi

Thu Aug 15 15:50:07 2024
+-----------------------------------------------------------------------------------------+
| NVIDIA-SMI 555.42.06              Driver Version: 555.42.06      CUDA Version: 12.6     |
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
|    2   N/A  N/A   1743472      C   ...r/.pixi/envs/default/bin/python3             0MiB |
+-----------------------------------------------------------------------------------------+ 
```

In addition to output file llama-cpp-server.out, the hostname and port information has also been saved in the `comment` field of the Slurm job and you can retrieve it from there, assuming your slurm job is 123456.

```
scontrol show job 123456 | sed -n 's/.*Comment=\([^[:space:]]*\).*/\1/p'

llama-cpp-server:cnode-12-2:33775
```

## Use openai api to connect to server

finally let's connect to our inference server form another system as a client. Right now this has to be run from a server that is a node  of the HPC system, for example the login node or the compute node we have an interactive session on. 

create a python script `llama-openai-client.py` with the content below 

```python
#!/usr/bin/env python3

import os, sys, openai

streaming=True

if len(sys.argv) < 3:
    print(f"Usage: {sys.argv[0]} <server:port> <prompt>")
    sys.exit(1)

server_port = sys.argv[1]
prompt = sys.argv[2]

client = openai.OpenAI(
    base_url=f"http://{server_port}/v1",
    api_key=os.getenv('OPENAI_API_KEY', os.getenv('USER') or os.getenv('USERNAME')),
)

def get_openai_response(client, prompt, streaming):
    try:
        response = client.chat.completions.create(
            messages=[
                {
                    "role": "user",
                    "content": prompt,
                }
            ],
            model="whatever",
            stream=streaming
        )
        # Process the stream
        if streaming:
            full_response = ""
            for chunk in response:
                chunk_message = chunk.choices[0].delta.content
                if chunk_message:
                    full_response += chunk_message
                    print(chunk_message, end='', flush=True)

            print()  # Print a newline after streaming is complete
            return full_response
        else:
            return response.choices[0].message.content
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return None

response = get_openai_response(client, prompt, streaming)
if not streaming: 
    print(response)
```

and then make it executable, and run it with the server and port name as the first argument and the prompt as the second command line argument. 

```
chmod +x llama-openai-client.py
llama-openai-client.py cnode-12-2:33775 "What is your home town?" 

  Meta, the company that created me, is headquartered in Menlo Park, California.

```

