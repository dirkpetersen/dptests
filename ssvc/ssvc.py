#! /usr/bin/env python3

import sys, os, argparse, subprocess, getpass
import yaml
import config_parser, ssh_utils, slurm_utils, load_balancer

# from config_parser import parse_configs 
# from ssh_utils import setup_ssh_key, SSHConnection
# from slurm_utils import submit_job, get_job_status
# from load_balancer import SimpleLoadBalancer

#ssh_utils.setup_ssh_key("exacloud.ohsu.edu") 
with ssh_utils.SSHConnection("exacloud.ohsu.edu") as ssh:
    ssh.setup_port_forwarding(8080, "exanode-11-34", 80)


def main():
    parser = argparse.ArgumentParser(description="Slurm Service Controller")
    parser.add_argument("config_file", help="Path to configuration YAML")
    args = parser.parse_args()

    configs = config_parser.parse_configs(args.config_file)

    # SSH setup
    if not ssh_utils.setup_ssh_key(configs['login_node']):
        print("SSH key setup failed. Exiting.")
        return 

    with ssh_utils.SSHConnection(configs['login_node']) as ssh:
        # Manage jobs, set up forwarding, etc.
        worker_configs = configs['workers']
        load_balancer = load_balancer.SimpleLoadBalancer(worker_configs)

        # ... (Job management, port forwarding logic) 

if __name__ == "__main__":
    main()


import ssh_utils 



