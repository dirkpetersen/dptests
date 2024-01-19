

document.querySelector('title').textContent = 'RFdiffusion: an open source method for structure generation';
**RFdiffusion: an open source method for structure generation** 


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



Rosetta Fold (RF) dissusion is an open source method for structure generation, with or without conditional information 
(a motif, target etc). It can perform motif scaffolding, unconditional protein generation, and other tasks.



### Reference:


* Joseph L. Watson, David Juergens, Nathaniel R. Bennett et al.  

*Broadly applicable and accurate protein design by integrating
structure prediction networks and diffusion generative models.*   

[arXiv preprint](https://www.biorxiv.org/content/10.1101/2022.12.09.519842v1.full), doi: https://doi.org/10.1101/2022.12.09.519842.

Documentation
	+ [RFdiffusion Github page](https://github.com/RosettaCommons/RFdiffusion)Important Notes
	+ Module Name: RFdissusion (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
	+ Unusual environment variables set
		- **RFDIFFUSION\_HOME**  installation directory
		- **RFDIFFUSION\_BIN**       executable directory
		- **RFDIFFUSION\_SRC**       source code directory
		- **RFDIFFUSION\_DATA**  sample data directory
Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=20g -c8 --gres=gpu:k80:1,lscratch:10**
[user@cn3335 ~]$ **module load RFdiffusion** 
[+] Loading singularity  3.10.5  on cn4338
[+] Loading RFdiffusion  1.1.0
[user@cn3335 ~]$ **git clone https://github.com/RosettaCommons/RFdiffusion**
[user@cn3335 ~]$ **cd RFdiffusion**
[user@cn3335 ~]$ **./scripts/run\_inference.py -h**
run_inference is powered by Hydra.

== Configuration groups ==
Compose your configuration from those groups (group=option)



== Config ==
Override anything in the config (foo.bar=value)

inference:
  input_pdb: null
  num_designs: 10
  design_startnum: 0
  ckpt_override_path: null
  symmetry: null
  recenter: true
  radius: 10.0
  model_only_neighbors: false
  output_prefix: samples/design
  write_trajectory: true
  scaffold_guided: false
  model_runner: SelfConditioning
  cautious: true
  align_motif: true
  symmetric_self_cond: true
  final_step: 1
  deterministic: false
  trb_save_ckpt_path: null
  schedule_directory_path: null
  model_directory_path: null
contigmap:
  contigs: null
  inpaint_seq: null
  provide_seq: null
  length: null
model:
  n_extra_block: 4
  n_main_block: 32
  n_ref_block: 4
  d_msa: 256
  d_msa_full: 64
  d_pair: 128
  d_templ: 64
  n_head_msa: 8
  n_head_pair: 4
  n_head_templ: 4
  d_hidden: 32
  d_hidden_templ: 32
  p_drop: 0.15
  SE3_param_full:
    num_layers: 1
    num_channels: 32
    num_degrees: 2
    n_heads: 4
    div: 4
    l0_in_features: 8
    l0_out_features: 8
    l1_in_features: 3
    l1_out_features: 2
    num_edge_features: 32
  SE3_param_topk:
    num_layers: 1
    num_channels: 32
    num_degrees: 2
    n_heads: 4
    div: 4
    l0_in_features: 64
    l0_out_features: 64
    l1_in_features: 3
    l1_out_features: 2
    num_edge_features: 64
  d_time_emb: null
  d_time_emb_proj: null
  freeze_track_motif: false
  use_motif_timestep: false
diffuser:
  T: 50
  b_0: 0.01
  b_T: 0.07
  schedule_type: linear
  so3_type: igso3
  crd_scale: 0.25
  partial_T: null
  so3_schedule_type: linear
  min_b: 1.5
  max_b: 2.5
  min_sigma: 0.02
  max_sigma: 1.5
denoiser:
  noise_scale_ca: 1
  final_noise_scale_ca: 1
  ca_noise_schedule_type: constant
  noise_scale_frame: 1
  final_noise_scale_frame: 1
  frame_noise_schedule_type: constant
ppi:
  hotspot_res: null
potentials:
  guiding_potentials: null
  guide_scale: 10
  guide_decay: constant
  olig_inter_all: null
  olig_intra_all: null
  olig_custom_contact: null
  substrate: null
contig_settings:
  ref_idx: null
  hal_idx: null
  idx_rf: null
  inpaint_seq_tensor: null
preprocess:
  sidechain_input: false
  motif_sidechain_input: true
  d_t1d: 22
  d_t2d: 44
  prob_self_cond: 0.0
  str_self_cond: false
  predict_previous: false
logging:
  inputs: false
scaffoldguided:
  scaffoldguided: false
  target_pdb: false
  target_path: null
  scaffold_list: null
  scaffold_dir: null
  sampled_insertion: 0
  sampled_N: 0
  sampled_C: 0
  ss_mask: 0
  systematic: false
  target_ss: null
  target_adj: null
  mask_loops: true
  contig_crop: null


Powered by Hydra (https://hydra.cc)
Use --hydra-help to view Hydra specific help

```

Download pretrained models:

```

[user@cn3335 ~]$ **bash scripts/download\_models.sh models**
...

```

Download sample data:

```

[user@cn3335 ~]$ **cp $RFDIFFUSION\_DATA/\* .**

```

If needed, edit the configuration file:

```

**config/inference/base.yaml**

```

Run RFdiffusion on the data, using the settings from the configuration file:

```

[user@cn3335 ~]$ **./scripts/run\_inference.py inference.output\_prefix=./ inference.input\_pdb=./sample.pdb 'contigmap.contigs=[10-40/a394-408/10-40]' +schedule\_directory\_path=./** 
[2023-06-05 10:23:41,241][__main__][INFO] - Found GPU with device_name Tesla K80. Will run RFdiffusion on Tesla K80
Reading models from /vf/users/user/RFdiffusion/RFdiffusion/rfdiffusion/inference/../../models
[2023-06-05 10:23:41,242][rfdiffusion.inference.model_runners][INFO] - Reading checkpoint from /vf/users/user/RFdiffusion/RFdiffusion/rfdiffusion/inference/../../models/Base_ckpt.pt
This is inf_conf.ckpt_path
/vf/users/user/RFdiffusion/RFdiffusion/rfdiffusion/inference/../../models/Base_ckpt.pt
Assembling -model, -diffuser and -preprocess configs from checkpoint
USING MODEL CONFIG: self._conf[model][n_extra_block] = 4
USING MODEL CONFIG: self._conf[model][n_main_block] = 32
USING MODEL CONFIG: self._conf[model][n_ref_block] = 4
USING MODEL CONFIG: self._conf[model][d_msa] = 256
USING MODEL CONFIG: self._conf[model][d_msa_full] = 64
USING MODEL CONFIG: self._conf[model][d_pair] = 128
USING MODEL CONFIG: self._conf[model][d_templ] = 64
USING MODEL CONFIG: self._conf[model][n_head_msa] = 8
USING MODEL CONFIG: self._conf[model][n_head_pair] = 4
USING MODEL CONFIG: self._conf[model][n_head_templ] = 4
USING MODEL CONFIG: self._conf[model][d_hidden] = 32
USING MODEL CONFIG: self._conf[model][d_hidden_templ] = 32
USING MODEL CONFIG: self._conf[model][p_drop] = 0.15
USING MODEL CONFIG: self._conf[model][SE3_param_full] = {'num_layers': 1, 'num_channels': 32, 'num_degrees': 2, 'n_heads': 4, 'div': 4, 'l0_in_features': 8, 'l0_out_features': 8, 'l1_in_features': 3, 'l1_out_features': 2, 'num_edge_features': 32}
USING MODEL CONFIG: self._conf[model][SE3_param_topk] = {'num_layers': 1, 'num_channels': 32, 'num_degrees': 2, 'n_heads': 4, 'div': 4, 'l0_in_features': 64, 'l0_out_features': 64, 'l1_in_features': 3, 'l1_out_features': 2, 'num_edge_features': 64}
USING MODEL CONFIG: self._conf[model][freeze_track_motif] = False
USING MODEL CONFIG: self._conf[model][use_motif_timestep] = True
USING MODEL CONFIG: self._conf[diffuser][T] = 50
USING MODEL CONFIG: self._conf[diffuser][b_0] = 0.01
USING MODEL CONFIG: self._conf[diffuser][b_T] = 0.07
USING MODEL CONFIG: self._conf[diffuser][schedule_type] = linear
USING MODEL CONFIG: self._conf[diffuser][so3_type] = igso3
USING MODEL CONFIG: self._conf[diffuser][crd_scale] = 0.25
USING MODEL CONFIG: self._conf[diffuser][so3_schedule_type] = linear
USING MODEL CONFIG: self._conf[diffuser][min_b] = 1.5
USING MODEL CONFIG: self._conf[diffuser][max_b] = 2.5
USING MODEL CONFIG: self._conf[diffuser][min_sigma] = 0.02
USING MODEL CONFIG: self._conf[diffuser][max_sigma] = 1.5
USING MODEL CONFIG: self._conf[preprocess][sidechain_input] = False
USING MODEL CONFIG: self._conf[preprocess][motif_sidechain_input] = True
USING MODEL CONFIG: self._conf[preprocess][d_t1d] = 22
USING MODEL CONFIG: self._conf[preprocess][d_t2d] = 44
USING MODEL CONFIG: self._conf[preprocess][prob_self_cond] = 0.5
USING MODEL CONFIG: self._conf[preprocess][str_self_cond] = True
USING MODEL CONFIG: self._conf[preprocess][predict_previous] = False
[2023-06-05 10:23:52,919][rfdiffusion.inference.model_runners][INFO] - Loading checkpoint.
[2023-06-05 10:23:58,119][rfdiffusion.diffusion][INFO] - Using cached IGSO3.
Successful diffuser __init__
[2023-06-05 10:23:58,199][__main__][INFO] - Making design ./_0
[2023-06-05 10:23:58,411][rfdiffusion.inference.model_runners][INFO] - Using contig: ['10-40/a394-408/10-40']
With this beta schedule (linear schedule, beta_0 = 0.04, beta_T = 0.28), alpha_bar_T = 0.00013696048699785024
[2023-06-05 10:23:58,462][rfdiffusion.inference.model_runners][INFO] - Sequence init: -----------------------LNETHFSDDIEQQAD-----------------------------------
[2023-06-05 10:24:04,076][rfdiffusion.inference.utils][INFO] - Sampled motif RMSD: 0.21
[2023-06-05 10:24:04,094][rfdiffusion.inference.model_runners][INFO] - Timestep 50, input to next step: -----------------------LNETHFSDDIEQQAD-----------------------------------
[2023-06-05 10:24:05,277][rfdiffusion.inference.utils][INFO] - Sampled motif RMSD: 0.19
[2023-06-05 10:24:05,281][rfdiffusion.inference.model_runners][INFO] - Timestep 49, input to next step: -----------------------LNETHFSDDIEQQAD-----------------------------------
[2023-06-05 10:24:06,866][rfdiffusion.inference.utils][INFO] - Sampled motif RMSD: 0.15
[2023-06-05 10:24:06,870][rfdiffusion.inference.model_runners][INFO] - Timestep 48, input to next step: -----------------------LNETHFSDDIEQQAD-----------------------------------
[2023-06-05 10:24:08,469][rfdiffusion.inference.utils][INFO] - Sampled motif RMSD: 0.14
[2023-06-05 10:24:08,473][rfdiffusion.inference.model_runners][INFO] - Timestep 47, input to next step: -----------------------LNETHFSDDIEQQAD-----------------------------------
[2023-06-05 10:24:10,060][rfdiffusion.inference.utils][INFO] - Sampled motif RMSD: 0.12
[2023-06-05 10:24:10,064][rfdiffusion.inference.model_runners][INFO] - Timestep 46, input to next step: -----------------------LNETHFSDDIEQQAD-----------------------------------
...
[2023-06-05 10:25:03,386][rfdiffusion.inference.utils][INFO] - Sampled motif RMSD: 0.18
[2023-06-05 10:25:03,390][rfdiffusion.inference.model_runners][INFO] - Timestep 2, input to next step: -----------------------LNETHFSDDIEQQAD-----------------------------------
[2023-06-05 10:25:05,602][__main__][INFO] - Finished design in 1.12 minutes
[2023-06-05 10:25:05,602][__main__][INFO] - Making design ./_1
[2023-06-05 10:25:05,664][rfdiffusion.inference.model_runners][INFO] - Using contig: ['10-40/a394-408/10-40']
With this beta schedule (linear schedule, beta_0 = 0.04, beta_T = 0.28), alpha_bar_T = 0.00013696048699785024
[2023-06-05 10:25:05,694][rfdiffusion.inference.model_runners][INFO] - Sequence init: ----------------------LNETHFSDDIEQQAD--------------------
[2023-06-05 10:25:06,601][rfdiffusion.inference.utils][INFO] - Sampled motif RMSD: 0.22
[2023-06-05 10:25:06,604][rfdiffusion.inference.model_runners][INFO] - Timestep 50, input to next step: ----------------------LNETHFSDDIEQQAD--------------------
...
[2023-06-05 10:25:51,862][__main__][INFO] - Finished design in 0.77 minutes
[2023-06-05 10:25:51,863][__main__][INFO] - Making design ./_2
[2023-06-05 10:25:51,925][rfdiffusion.inference.model_runners][INFO] - Using contig: ['10-40/a394-408/10-40']
With this beta schedule (linear schedule, beta_0 = 0.04, beta_T = 0.28), alpha_bar_T = 0.00013696048699785024
[2023-06-05 10:25:51,955][rfdiffusion.inference.model_runners][INFO] - Sequence init: -------------------------------LNETHFSDDIEQQAD-------------
[2023-06-05 10:25:52,886][rfdiffusion.inference.utils][INFO] - Sampled motif RMSD: 0.20
[2023-06-05 10:25:52,889][rfdiffusion.inference.model_runners][INFO] - Timestep 50, input to next step: -------------------------------LNETHFSDDIEQQAD-------------
...
[2023-06-05 10:32:34,898][rfdiffusion.inference.utils][INFO] - Sampled motif RMSD: 0.13
[2023-06-05 10:32:34,902][rfdiffusion.inference.model_runners][INFO] - Timestep 2, input to next step: -------------------LNETHFSDDIEQQAD------------------------
[2023-06-05 10:32:36,721][__main__][INFO] - Finished design in 0.79 minutes

```

End the interactive session:

```

[user@cn3335 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```
