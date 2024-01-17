
vi .local/lib/python3.11/site-packages/easybuild/easyblocks/generic/conda.py

fix 

#        cmd = "%s config --add create_default_packages setuptools" % conda_cmd
#        run_cmd(cmd, log_all=True, simple=True)


and 


            #cmd = "%s %s create --force -y -p %s %s"
            cmd = "%s %s create -y -p %s %s" % (self.cfg['preinstallopts'], conda_cmd,
