

document.querySelector('title').textContent = "TEMPLATE";
CHAP on Biowulf


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



 A new tool for the functional annotation of novel ion channel structures that provides information on the biophysical properties of the ion permeation pathway by utilising molecular dynamics simulations.



### References:


* Gianni Klesse, Shanlin Rao, Mark S.P. Sansom, Stephen J. Tucker
 [**CHAP: A Versatile Tool for the Structural and Functional Annotation of Ion Channel Pores**](https://doi.org/10.1016/j.jmb.2019.06.003)
*Volume 431, Issue 17, 9 August 2019, Pages 3353-3365*


Documentation
* [Chap Main Site](http://www.channotation.org/)


Important Notes
* Module Name: chap (see [the modules page](/apps/modules.html) for more information)
	+ TEMPLATE\_HOME* Example files in /usr/local/apps/chap/chap-version\_0\_9\_1/examples



Getting Started
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
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

[user@cn3144 ~]$ **module load chap**
[+] Loading chap  0.9.1  on cn3144
[+] Loading singularity  3.10.5  on cn3144

[user@cn3144 ~]$**chap -h**

CCCCCC  HH     HH    AAA    PPPPPPPP
CC    CC HH     HH   AA AA   PP     PP
CC       HH     HH  AA   AA  PP     PP
CC       HHHHHHHHH AA     AA PPPPPPPP
CC       HH     HH AAAAAAAAA PP
CC    CC HH     HH AA     AA PP
 CCCCCC  HH     HH AA     AA PP

The Channel Annotation Package, version 0.9.1

SYNOPSIS

chap [-f [<.xtc/.trr/...>]] [-s [<.tpr/.gro/...>]] [-n [<.ndx>]] [-b ]
[-e ] [-dt ] [-tu ] [-fgroup ]
[-xvg ] [-sf ] [-selrpos ] [-seltype ]
[-sel-pathway ] [-sel-solvent ]
[-out-filename ] [-out-num-points ]
[-out-extrap-dist ] [-out-grid-dist ]
[-out-vis-tweak ] [-[no]out-detailed] [-pf-method ]
[-pf-vdwr-database ] [-pf-vdwr-fallback ]
[-pf-vdwr-json ] [-pf-align-method ]
[-pf-probe-step ] [-pf-max-free-dist ]
[-pf-max-probe-steps ] [-pf-sel-ipp ]
[-pf-init-probe-pos ] [-pf-chan-dir-vec ]
[-pf-cutoff ] [-sa-seed ] [-sa-max-iter ]
[-sa-init-temp ] [-sa-cooling-fac ] [-sa-step ]
[-nm-max-iter ] [-nm-init-shift ] [-pm-pl-margin ]
[-pm-pf-sel ] [-de-method ] [-de-res ]
[-de-bandwidth ] [-de-bw-scale ] [-de-eval-cutoff ]
[-hydrophob-database ] [-hydrophob-fallback ]
[-hydrophob-json ] [-hydrophob-bandwidth ]

DESCRIPTION

CHAP finds pores in biological macromolecules like ion channels and determines
the hydration state of these permeation pathways. It can operate on both
individual structures and on molecular dynamics trajectories. Visit
www.channotation.org for the full documentation.

OPTIONS

Options to specify input files:

-f [<.xtc/.trr/...>] (traj.xtc) (Opt.)
Input trajectory or single configuration: xtc trr cpt gro g96 pdb
tng
-s [<.tpr/.gro/...>] (topol.tpr) (Opt.)
Input structure: tpr gro g96 pdb brk ent
-n [<.ndx>] (index.ndx) (Opt.)
Extra index groups

Other options:

-b  (0)
First frame (ps) to read from trajectory
-e  (0)
Last frame (ps) to read from trajectory
-dt  (0)
Only use frame if t MOD dt == first time (ps)
-tu  (ps)
Unit for time values: fs, ps, ns, us, ms, s
-fgroup 
Atoms stored in the trajectory file (if not set, assume first N
atoms)
-xvg  (xmgrace)
Plot formatting: none, xmgrace, xmgr
-sf 
Provide selections from files
-selrpos  (atom)
Selection reference positions: atom, res\_com, res\_cog, mol\_com,
mol\_cog, whole\_res\_com, whole\_res\_cog, whole\_mol\_com,
whole\_mol\_cog, part\_res\_com, part\_res\_cog, part\_mol\_com,
part\_mol\_cog, dyn\_res\_com, dyn\_res\_cog, dyn\_mol\_com, dyn\_mol\_cog
-seltype  (atom)
Default selection output positions: atom, res\_com, res\_cog,
mol\_com, mol\_cog, whole\_res\_com, whole\_res\_cog, whole\_mol\_com,
whole\_mol\_cog, part\_res\_com, part\_res\_cog, part\_mol\_com,
part\_mol\_cog, dyn\_res\_com, dyn\_res\_cog, dyn\_mol\_com, dyn\_mol\_cog
-sel-pathway 
Reference group that defines the permeation pathway (usually
'Protein')
-sel-solvent 
Group of small particles to calculate density of (usually 'Water')
-out-filename  (output)
File name for output files without file extension.
-out-num-points  (1000)
Number of spatial sample points that are written to the JSON output
file.
-out-extrap-dist  (0)
Extrapolation distance beyond the pathway endpoints for both JSON
and OBJ output.
-out-grid-dist  (0.15)
Controls the sampling distance of vertices on the pathway surface
which are subsequently interpolated to yield a smooth surface. Very
small values may yield visual artifacts.
-out-vis-tweak  (0.1)
Visual tweaking factor that controls the smoothness of the pathway
surface in the OBJ output. Varies between -1 and 1 (exculisvely),
where larger values result in a smoother surface. Negative values
may result in visualisation artifacts.
-[no]out-detailed (no)
If true, CHAP will write detailed per-frame information to a
newline delimited JSON file including original probe positions and
spline parameters. This is mostly useful for debugging.
-pf-method  (inplane\_optim)
Path finding method. The default inplane\_optim implements the
algorithm used in the HOLE programme, where the position of a probe
sphere is optimised in subsequent parallel planes so as to maximise
its radius. The alternative naive\_cylindrical simply uses a
cylindrical volume as permeation pathway.: cylindrical,
inplane\_optim
-pf-vdwr-database  (hole\_simple)
Database of van-der-Waals radii to be used in pore finding:
hole\_amberuni, hole\_bondi, hole\_hardcore, hole\_simple, hole\_xplor,
user
-pf-vdwr-fallback  (nan)
Fallback van-der-Waals radius for atoms that are not listed in
van-der-Waals radius database. Unless this is set to a positive
value, an error will be thrown if a pathway-forming atom has no
associated van-der-Waals radius in the database.
-pf-vdwr-json 
JSON file with user defined van-der-Waals radii. Will be ignored
unless -pf-vdwr-database is set to 'user'.
-pf-align-method  (ipp)
Method for aligning pathway coordinates across time steps: none,
ipp
-pf-probe-step  (0.1)
Step length for probe movement.
-pf-max-free-dist  (1)
Maximum radius of pore.
-pf-max-probe-steps  (10000)
Maximum number of steps the probe is moved in either direction.
-pf-sel-ipp 
Selection of atoms whose COM will be used as initial probe
position. If not set, the selection specified with 'sel-pathway'
will be used.
-pf-init-probe-pos  (nan nan nan)
Initial position of probe in probe-based pore finding algorithms.
If set explicitly, it will overwrite the COM-based initial position
set with the ippSelflag.
-pf-chan-dir-vec  (0 0 1)
Channel direction vector. Will be normalised to unit vector
internally.
-pf-cutoff  (nan)
Cutoff for distance searches in path finding algorithm. A value of
zero or less means no cutoff is applied. If unset, an appropriate
cutoff is determined automatically.
-sa-seed  (0)
Seed used in pseudo random number generation for simulated
annealing. If not set explicitly, a random seed is used.
-sa-max-iter  (0)
Number of cooling iterations in one simulated annealing run.
-sa-init-temp  (0.1)
Simulated annealing initial temperature.
-sa-cooling-fac  (0.98)
Simulated annealing cooling factor.
-sa-step  (0.001)
Step length factor used in candidate generation.
-nm-max-iter  (100)
Number of Nelder-Mead simplex iterations in path finding algorithm.
-nm-init-shift  (0.1)
Distance of vertices in initial Nelder-Mead simplex.
-pm-pl-margin  (0.75)
Margin for determining pathway lining residues. A residue is
considered to be pathway lining if it is no further than the local
path radius plus this margin from the pathway's centre line.
-pm-pf-sel  (name CA)
Selection string that determines the group of atoms in each residue
whose centre of geometry's distance from the centre line is used to
determine if a residue is pore-facing.
-de-method  (kernel)
Method used for estimating the probability density of the solvent
particles along the permeation pathway: histogram, kernel
-de-res  (0.01)
Spatial resolution of the density estimator. In case of a
histogram, this is the bin width, in case of a kernel density
estimator, this is the spacing of the evaluation points.
-de-bandwidth  (-1)
Bandwidth for the kernel density estimator. Ignored for other
methods. If negative or zero, bandwidth will be determined
automatically to minimise the asymptotic mean integrated squared
error (AMISE).
-de-bw-scale  (1)
Scaling factor for the band width. Useful to set a bandwidth
relative to the AMISE-optimal value.
-de-eval-cutoff  (5)
Evaluation range cutoff for kernel density estimator in multiples
of bandwidth. Ignored for other methods. Ensures that the density
falls off smoothly to zero outside the data range.
-hydrophob-database  (wimley\_white\_1996)
Database of hydrophobicity scale for pore forming residues:
hessa\_2005, kyte\_doolittle\_1982, monera\_1995, moon\_2011,
wimley\_white\_1996, zhu\_2016, memprotmd, user
-hydrophob-fallback  (nan)
Fallback hydrophobicity for residues in the pathway defining group.
If unset (nan), residues missing in the database will cause an
error.
-hydrophob-json 
JSON file with user defined hydrophobicity scale. Will be ignored
unless -hydrophobicity-database is set to 'user'.
-hydrophob-bandwidth  (0.35)
Bandwidth for hydrophobicity kernel.


Thank you for using CHAP - The Channel Annotation Package!

[user@cn3144 ~]$ 

```


Example
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
To compute a radius profile for this channel structure



```

[user@cn3144 ~]$**cp -r /usr/local/apps/chap/chap-version\_0\_9\_1/examples/example-01 .** 
[user@cn3144 ~]$ **cd example-01/**
[user@cn4282 example-01]$**chap -f 4pirtm.pdb -s 4pirtm.pdb**

CCCCCC  HH     HH    AAA    PPPPPPPP
CC    CC HH     HH   AA AA   PP     PP
CC       HH     HH  AA   AA  PP     PP
CC       HHHHHHHHH AA     AA PPPPPPPP
CC       HH     HH AAAAAAAAA PP
CC    CC HH     HH AA     AA PP
 CCCCCC  HH     HH AA     AA PP

The Channel Annotation Package, version 0.9.1


WARNING: Masses and atomic (Van der Waals) radii will be guessed
based on residue and atom names, since they could not be
definitively assigned from the information in your input
files. These guessed numbers might deviate from the mass
and radius of the atom type. Please check the output
files if necessary.

Available static index groups:
Group  0 "System" (175315 atoms)
Group  1 "Protein" (12113 atoms)
Group  2 "Protein-H" (5838 atoms)
Group  3 "C-alpha" (722 atoms)
Group  4 "Backbone" (2166 atoms)
Group  5 "MainChain" (2878 atoms)
Group  6 "MainChain+Cb" (3575 atoms)
Group  7 "MainChain+H" (3570 atoms)
Group  8 "SideChain" (8543 atoms)
Group  9 "SideChain-H" (2960 atoms)
Group 10 "Prot-Masses" (12113 atoms)
Group 11 "non-Protein" (163202 atoms)
Group 12 "Other" (61774 atoms)
Group 13 "POPC" (61774 atoms)
Group 14 "NA" (163 atoms)
Group 15 "CL" (168 atoms)
Group 16 "Water" (101097 atoms)
Group 17 "SOL" (101097 atoms)
Group 18 "non-Water" (74218 atoms)
Group 19 "Ion" (331 atoms)
Group 20 "POPC" (61774 atoms)
Group 21 "NA" (163 atoms)
Group 22 "CL" (168 atoms)
Group 23 "Water_and_ions" (101428 atoms)
Specify a selection for option 'sel-pathway'
(Reference group that defines the permeation pathway (usually 'Protein') ):
(one per line,  for status/groups, 'help' for help)
> 1
Selection '1' parsed

Reading frame 0 time 0.000 '', 175315 atoms
Last frame 0 time 0.000
Analyzed 1 frames, last time 0.000

Forming time averages, 100% complete


Thank you for using CHAP - The Channel Annotation Package!


```

For more examples please see the [Examples page](http://www.channotation.org/docs/annotation_example/) 








