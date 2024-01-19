

document.querySelector('title').textContent = 'Mathematica on Biowulf ';
Mathematica on Biowulf 


|  |
| --- |
| 
Quick Links
[Run Mathematica GUI Interactively](#gui)
[Run Mathematica without GUI](#nogui)
[Run Mathematica using math file](#math)
[Run Mathematica jobs on Biowulf](#jobs)
[Troubleshooting](#trouble)
[Documentation](#docs)
 |


Mathematica is a fully integrated environment for
technical and scientific computing. Mathematica combines numerical and
symbolic computation, visualization, and programming in a single,
flexible interactive system.
The Mathematica system is very broad, and provides a systematic
interface to all sorts of computations, from traditional numeric and
symbolic computation, to visualization, to data format conversion, and
the creation of user interfaces.


Running the Mathematica GUI Interactively
To run the Mathematica graphics interface, an X-Windows connection is required (preferably using [Nomachine](/docs/nx.html)). Open an X-Windows connection to biowulf.nih.gov, and check that it is working by typing 'xclock' at the prompt.



Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job
[user@cn3144 ~]$ **module load mathematica**
[+] Loading mathematica 12.1  on cn3144 

[user@cn3144 ~]$ **mathematica &**

```

You should see the Mathematica Welcome Screen as below:


![Mathematica Welcome Screen](/images/mathematica_12.1_splash_screen.png)


The main Mathematica Notebook window will appear when you open a file by clicking Open or a file listed on the Welcome screen.


To run a demo, in the main Mathematica window, click Open. Enter the file name /usr/local/apps/Mathematica/12.1/Documentation/English/System/HowTos/CustomizeMyGraphics.nb. You should see a window appear with a mathematica notebook example:


![Mathematica demo image](/images/mathematica_12.1_example.png)


Run Mathematica interactively without the GUI
It is also possible to run Mathematica interactively on the command-line, without the GUI. This would be useful if you do not wish to use Xwindows, and is obviously most useful for calculations rather than interactive graphics. Sample session (user input in bold):



Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job
[user@cn3144 ~]$ **module load mathematica**
[+] Loading mathematica 12.1  on cn3144 

[user@cn3144 ~]$ **math**
Mathematica 12.1.0 Kernel for Linux x86 (64-bit)
Copyright 1988-2020 Wolfram Research, Inc.

In[1]:=  **data=ReadList["rmsdplot.list",{Number,Number,Number}];**

In[2]:= [..enter other Mathematica commands..]

In[22]:= **Exit[]**
[user@cn3144 ~]$
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

Run Mathematica using a math command file
You can insert Mathematica commands into a math command file, and run
that file interactively at the prompt. This is convenient if you
perform the same Mathematica tasks frequently. Let us use the following 
Math command file test.m.



```

------------------------ test.m -------------------------------
A = Sum[i, {i,1,100}]
B = Mean[{25, 36, 22, 16, 8, 42}]
Answer = A + B
Print["Answer=", Answer]
Exit[]
-------------------- test.m -----------------------------------

```

  

Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):

```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job
[user@cn3144 ~]$ **module load mathematica**
[+] Loading mathematica 12.1  on cn3144 

[user@cn3144 ~]$ **math**
Mathematica 12.1.0 Kernel for Linux x86 (64-bit)
Copyright 1988-2020 Wolfram Research, Inc.

In[1]:= **<<test.m**

       30449
Answer=-----
         6

[user@cn3144 ~]$

```

Alternatively, you can type directly at the command prompt:


```

[user@cn3144~] **math -noprompt -run "<<test.m"**
Answer=30449/6
[user@cn3144 ~]

```

Running Mathematica Jobs on Biowulf
Mathematica provides a seamless, integrated and automated environment for parallel computing. To submit a Mathematica job to the batch system, you must create a .m file containing the mathematica commands. The .m script will be passed to the math command in an sbatch file. 

The above example can be sent to batch with the following script:


```

------------------------ test.sh-------------------------------
#!/bin/bash

module load mathematica

math -run < test.m > outputfile.out
--------------------end test.sh -------------------------------

```

Submit with **sbatch test.sh**
Troubleshooting
If Mathematica prompts you for an Activation Key:

![Mathematica Activation Screen](/images/mathematica_activation_screenshot.jpg)



Click on "Other ways to activate", and set the licenses server to **badmin.cit.nih.gov**.


Documentation
[Wolfram Mathematica Documentation Center](https://reference.wolfram.com/language)
  
[Mathematica FAQs](https://www.wolfram.com/FAQs/)














