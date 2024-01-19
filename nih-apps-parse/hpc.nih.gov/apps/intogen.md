

document.querySelector('title').textContent = 'intogen: dentifies cancer genes and pinpoints their putative mechanism of action across tumor types';
**intogen: dentifies cancer genes and pinpoints their putative mechanism of action across tumor types**


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



- IntOGen is a framework for automatic and comprehensive knowledge extraction based on mutational data from sequenced tumor samples from patients. The framework identifies cancer genes and pinpoints their putative mechanism of action across tumor types.

- Given a dataset of somatic point mutations from a cohort of tumor samples, intOGen first pre-processes the input mutations, next it runs seven different methods for cancer driver gene identification and, finally, it combines the output of these methods to produce a compendium of driver genes and a repository of the mutational features that can be used to explain their mechanisms of action.

- They have manually downloaded and annotated tumor samples from different sources. Specifically, we have used cBioPortal, pediatric cBioPortal, ICGC, TCGA, PCAWG, Hartwig Medical Foundation, TARGET, St. Jude and literature gathered sequencing projects projects.

### References:



Francisco Martínez-Jiménez, Ferran Muiños, Inés Sentís, Jordi Deu-Pons, Iker Reyes-Salazar, Claudia Arnedo-Pac, Loris Mularoni, Oriol Pich, Jose Bonet, Hanna Kranas, Abel Gonzalez-Perez & Nuria Lopez-Bigas 
  

*A compendium of mutational cancer driver gene*   

[PubMed](https://pubmed.ncbi.nlm.nih.gov/32778778/)
[Nature Reviews Caner](https://www.nature.com/articles/s41568-020-0290-x) 2020



Documentation
* [intogen bitbucket Page](https://bitbucket.org/intogen/intogen-plus/src/master/)
* [intogen doc](https://intogen.readthedocs.io/en/latest/index.html)


Important Notes
* Module Name: intogen (see [the modules page](/apps/modules.html) for more information)
* Unusual environment variables set
	+ **INTOGEN\_TEST\_DATA** sample data for running intogen



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive -c 14 --mem=20g --gres=lscratch:10** 
[user@cn3144 ~]$ **module load intogen** 
[+] Loading singularity  3.10.5  on cn4304
[+] Loading java 17.0.2  ...
[+] Loading nextflow 22.10.2
[+] Loading intogen  2023.1

```


```

[user@cn3144 ]$ **cp -r $INTOGEN\_TEST\_DATA/\* .** 
[user@cn3144 ]$ **cp -r $INTOGEN\_CONF/\* .**
[user@cn3144 ]$ **nextflow run intogen.nf -profile local --input test/ --output ./output**
N E X T F L O W  ~  version 22.10.4
Launching `intogen.nf` [curious_ride] DSL1 - revision: eb5992f5c0
executor >  local (33)
[3a/3b8d38] process > ParseInput (Parse input test)                                                       [100%] 1 of 1 ✔
[4e/dfea6a] process > LoadCancer (Load cancer type CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)                       [100%] 1 of 1 ✔
[9d/ba5144] process > LoadPlatform (Load sequencing platform CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)             [100%] 1 of 1 ✔
[81/c96e65] process > LoadGenome (Load reference genome CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)                  [100%] 1 of 1 ✔
[7b/fd10a2] process > ProcessVariants (Process variants CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)                  [100%] 1 of 1 ✔
[de/4a398e] process > FormatSignature (Prepare for signatures CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)            [100%] 1 of 1 ✔
[d5/b8a0e0] process > Signature (Signatures CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)                              [100%] 1 of 1 ✔
[2a/da6345] process > FormatFML (Prepare for FML CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)                         [100%] 1 of 1 ✔
[42/77362b] process > OncodriveFML (OncodriveFML CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)                         [100%] 1 of 1 ✔
[4c/143f72] process > FormatCLUSTL (Prepare for CLUSTL CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)                   [100%] 1 of 1 ✔
[8e/9258f2] process > OncodriveCLUSTL (OncodriveCLUSTL CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)                   [100%] 1 of 1 ✔
[11/038f2f] process > FormatDNDSCV (Prepare for DNDSCV CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)                   [100%] 1 of 1 ✔
[11/79a963] process > dNdScv (dNdScv CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)                                     [100%] 1 of 1 ✔
[5c/e60332] process > FormatVEP (Prepare for VEP CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)                         [100%] 1 of 1 ✔
[40/34e5f0] process > VEP (VEP CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)                                           [100%] 1 of 1 ✔
[98/cdc7c0] process > ProcessVEPoutput (Process vep output CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)               [100%] 1 of 1 ✔
[10/43e0de] process > FilterNonSynonymous (Filter non synonymus CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)          [100%] 1 of 1 ✔
[b2/3725ce] process > FormatSMRegions (Prepare for SMRegions CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)             [100%] 1 of 1 ✔
[21/c1f73f] process > SMRegions (SMRegions CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)                               [100%] 1 of 1 ✔
[d2/2a5cf1] process > FormatCBaSE (Prepare for CBaSE CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)                     [100%] 1 of 1 ✔
[70/232727] process > CBaSE (CBaSE CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)                                       [100%] 1 of 1 ✔
[2e/5fec7b] process > FormatMutPanning (Prepare for MutPanning CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)           [100%] 1 of 1 ✔
[4e/7fbafe] process > MutPanning (MutPanning CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)                             [100%] 1 of 1 ✔
[87/a65c42] process > FormatHotMAPS (Prepare for HotMAPS CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)                 [100%] 1 of 1 ✔
[a8/96f071] process > HotMAPS (HotMAPS CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)                                   [100%] 1 of 1 ✔
[4d/1c88d2] process > Combination (Combination CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)                           [100%] 1 of 1 ✔
[da/7d1615] process > FormatdeconstructSigs (Prepare for deconstructSigs CBIOP_WXS_CBIOPORTAL_PRAD_BROAD) [100%] 1 of 1 ✔
[48/a23b2d] process > deconstructSigs (deconstructSigs CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)                   [100%] 1 of 1 ✔
[8f/07cd96] process > CohortCounts (Count variants CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)                       [100%] 1 of 1 ✔
[d5/c8d224] process > CohortSummary (Count variants null)                                                 [100%] 1 of 1 ✔
[63/0f19ec] process > MutationsSummary (Mutations)                                                        [100%] 1 of 1 ✔
[f6/aebe60] process > DriverDiscovery (Driver discovery CBIOP_WXS_CBIOPORTAL_PRAD_BROAD)                  [100%] 1 of 1 ✔
[69/a1269d] process > DriverSummary (Driver summary)                                                      [100%] 1 of 1 ✔
Completed at: 22-Feb-2023 15:53:34
Duration    : 3h 36m 50s
CPU hours   : 15.6
Succeeded   : 33


```







