import ssh_utils 

#ssh_utils.setup_ssh_key("exacloud.ohsu.edu") 

with ssh_utils.SSHConnection("exacloud.ohsu.edu") as ssh:
    ssh.setup_port_forwarding(8080, "exanode-11-34", 80)
