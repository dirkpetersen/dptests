# Slurm Service 

Feed this info to chat gpt or other software making things: 

I want to run a pool of rest apis on one or multiple slurm compute nodes as a backend service to another rest api which runs on an application server outside the slurm cluster. 

The slurm compute node can only be reached by the slurm login node but not by the application server and the only way for the application server to reach the compute node is ssh port forwarding using the login node as a conduit. 

The configuration details of the rest services on the compute nodes (path, command line, port, working directory, # cpu, memory GiB, # gpu and type (optional), local disk space GiB, and shared disk space GiB requirements) should be configured in a service yaml file that is located at a url (either https:// or file:///). 

Write a python tool with a cli named ssvc <config.yaml> . This config yaml file is different from the service yaml file in that it describes the slurm infrastructure (login node fqdn, slurm partition, slurm qos, gpu request option (e.g. --gres), local disk request option (e.g. --gres) options for getting gpu types and local scratch spaces, number of REST processes to run and the location of the service yaml file).

When ssvc <config.yaml> is executed it will connect to the login node via ssh (please use plain ssh and not paramiko), check how many of the desired workers are activelty running and request the number of workers that are missing. if the number of workers is complete implement port forwarding so that each of the workers can be reached on a different port on 127.0.0.1 on the application server. The rest api on the application server should implement load balancing as well as failover using all the worker nodes 

After the rest service on the compute node has been started using sbatch it should run a command like this so other applications can find where these services are running by using the squeue command and parsing the output of the comment 

f'scontrol update job {SLURM_JOB_ID} comment="app={appname}|job={jobname}|host={hostname}|port={portnum}"'

please check this url to see how this could be implemented 
https://gist.github.com/Delaunay/8c866a81cd696ca4cc01df26d6849764

The ssvc tool should connect with the same username it runs on the applocation server to the loginnode. if it does not find an ssh key that has the name ssvc in the file name it should automatically generate a passwordless ssh key and scp that public key to the login node and prompt for a password and then add that public key to authorized keys on the login node 

General coding instructions: 
- always use shebang #! /usr/bin/env python3
- needs to be compatible with Python 3.9 and newer 
- put all import statements into as few lines as possible, line length not to exceed 80 chars 
- don't mix and match import statements from the standard library with external modules
- generate a new requirements.txt whenever you propose a new external module 
- use f-strings with single quotes whenever possible 
- use async io whenever feasible 
- use os.path.expanduser to get a home directory 
- use print statements with flush so they 
- ensure that python subprocess will always print to STDOUT/STDERR right away 

SSH instructions 
- use SSH_OPTIONS = "-o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no -o LogLevel=ERROR"
- never use paramiko, only subprocess to standard openssh 



Implement all the code for this python project 
