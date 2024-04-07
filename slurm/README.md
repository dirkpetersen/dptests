
Write me a tool 'tsqueue' with a TUI in python that lists the output of squeue using https://github.com/Textualize/textual. Using the DataTable widget is one option 

 scope: - all jobs running under my default slurm account, my own jobs first 
        - if I don have a default slurm account just list my own jobs 
        - read the columns from env var export SQUEUE_FORMAT
        - by default SQUEUE_FORMAT should be "%.18i %.4P %.12j %.8u %.2t %.10M %.10L %.3D %.3C %.9b %.4m %R"
        - if a job is not running show a lay explanation for it. 

I need to navigate the job-list with up down keys, when hitting enter on a job it should show all possible details that this job might have available, here is some sample code that retrieves those details, it is bad code though 

https://github.com/dirkpetersen/hpctoys/blob/main/bin/s-jobinfo

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


