
vi .local/lib/python3.11/site-packages/easybuild/easyblocks/generic/conda.py

fix 

#        cmd = "%s config --add create_default_packages setuptools" % conda_cmd
#        run_cmd(cmd, log_all=True, simple=True)


and 


            #cmd = "%s %s create --force -y -p %s %s"
            cmd = "%s %s create -y -p %s %s" % (self.cfg['preinstallopts'], conda_cmd,


thios works :

micromamba create -n test319 python=3.8 -c conda-forge


diff -ruN dummy micromamba-wrapper.sh > micromamba-wrapper.patch

patch < micromamba-wrapper.patch

change :

--- micromamba  1969-12-31 16:00:00.000000000 -0800
+++ micromamba  2024-01-17 10:41:49.921401842 -0800
@@ -0,0 +1,51 @@
