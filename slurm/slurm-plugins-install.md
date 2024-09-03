# Installing several slurm Plugins


First let's create a separate plugin dir for custom 3rd party plugins 

```
mkdir -p /etc/slurm/plugins/23.11.10
```

and add this dir to the PluginDir setting in slurm.conf. Since spank plugins need to be recompiled for each version of slurm it is recommended to use a dedicated folder per Slurm version 

```
PluginDir=/usr/lib64/slurm:/etc/slurm/plugins/23.11.10
#PlugStackConfig=
```


## auto_tmpdir 

Please follow these steps to build auto_tmpdir

```
git clone git@github.com:University-of-Delaware-IT-RCI/auto_tmpdir.git
cd auto_tmpdir
mkdir build 
cd build 

cmake -DCMAKE_BUILD_TYPE=Release \
  -DSLURM_MODULES_DIR=/usr/lib64/slurm \  
  ..

make 

cp auto_tmpdir.so /etc/slurm/plugins/23.11.10

mkdir /var/tmp/auto_tmpdir_cache
chown slurm /var/tmp/auto_tmpdir_cache


```

and finally we edit /etc/slurm/plugstack.conf and add a line 

```
optional    auto_tmpdir.so  mount=/tmp mount=/var/tmp local_prefix=/mnt/scratch/tmpdir- state_dir=/var/tmp/auto_tmpdir_cache
```

and restart slurmd 

```
systemctl restart slurmd 
```

and test as rocky user and see the empty tmpdirs 

```
su - rocky 
srun --pty bash 

ls $TMPDIR /tmp /var/tmp /dev/shm
/dev/shm:

/tmp:

/tmp:

/var/tmp:
```

## auto_tmpdir with shared tmpdir

This does not seem to work : 

cmake -DCMAKE_BUILD_TYPE=Release \
   -DAUTO_TMPDIR_ENABLE_SHARED_TMPDIR=On \
   -DAUTO_TMPDIR_DEFAULT_SHARED_PREFIX=/arc/scratch1/jobs \
   -DSLURM_MODULES_DIR=/usr/lib64/slurm \
   ..


in /etc/slurm/plugstack.conf 

```
optional    auto_tmpdir.so  mount=/tmp mount=/var/tmp local_prefix=/mnt/scratch/tmpdir- shared_prefix=/arc/scratch1/slurm/job- no_rm_shared_only state_dir=/var/tmp/auto_tmpdir_cache
```
