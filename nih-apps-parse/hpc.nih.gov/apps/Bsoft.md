

document.querySelector('title').textContent = ' Bsoft on Biowulf ';

Bsoft on Biowulf 



|  |
| --- |
| 
Quick Links
[Batch job on Biowulf](#batch)
[Swarm of jobs](#swarm)
[Interactive job on Biowulf](#int)
[Documentation](#doc)
 |


[![](/images/Bsoft.jpg)](http://bsoft.ws)Bsoft is a collection of programs and a platform for development of software for image and molecular processing in structural biology. Problems in structural biology are approached with a highly modular design, allowing fast development of new algorithms without the burden of issues such as file I/O. It provides an easily accessible interface, a resource that can be and has been used in other packages.

Bsoft was developed by groups at the University of Basel, NIAMS, and Caltech.   

[Bsoft website](http://bsoft.ws)


Batch job on Biowulf

The Bsoft programs have standard Unix command-line syntax, and so are ideally suited to running as a batch job or as a swarm of batch jobs. 

Sample Bsoft batch script:

```

#!/bin/bash

module load Bsoft

cd /data/user/myjob
bimg -verbose 7 -truncate 0,100 input.img output.img

```


This job would be submitted with

```

sbatch  jobscript

```



This job would be submitted to a single core (2 hyperthreaded CPUs) and 4 GB of memory. If more than the default 4 GB of 
memory is required, you would submit with:

```

sbatch --mem=#g  jobscript

```

where '#' is the number of GigaBytes of memory required. 

Swarm of jobs

Biowulf users will typically want to process large numbers of image files via Bsoft. This is most easily done using the swarm utility. 


Now set up a swarm command file along the following lines:

```

# this file is called swarmfile

bimg -verbose 7 -truncate 0,100 input1.img output1.img
bimg -verbose 7 -truncate 0,100 input2.img output2.img
bimg -verbose 7 -truncate 0,100 input3.img output3.img
bimg -verbose 7 -truncate 0,100 input4.img output4.img
bimg -verbose 7 -truncate 0,100 input5.img output5.img
[...]

```


Submit this with

```

biowulf% swarm -f swarmfile --module Bsoft/2.0.0

```


If each Bsoft process in the swarm command file requires more than 4 GB of memory, use

```

biowulf% swarm -g # -f swarmfile --module Bsoft/1.9.0

```

where # is the number of Gigabytes required by each process.


Interactive job

To run Bsoft interactively on Biowulf, allocate an interactive session and run the process there.

Sample session:

```

biowulf% sinteractive

cn0203% module load bsoft

cn0203% ...bsoft commands ...

cn0203% exit

```

Documentation

[Bsoft website](http://lsbr.niams.nih.gov/Bsoft/bsoft.html)






































