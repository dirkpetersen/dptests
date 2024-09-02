#! /usr/bin/python3

import subprocess, os, re, grp, click

class BigBadClass:
    def __init__(self):
        self.username = os.getenv("USER")
        self.user_groups = [g.gr_name for g in grp.getgrall() if self.username in g.gr_mem]
        print(f"User: {self.username}, Groups: {self.user_groups}")

        self.partitions_info = self._parse_scontrol_output(self._run_command('scontrol --oneline show partition'))
        print(f"Partitions Info: {self.partitions_info}")

        self.sinfo_data = self._parse_sinfo_output(self._run_command('sinfo -o "%P|%l|%C"'))
        print(f"Sinfo Data: {self.sinfo_data}")

        self.user_permissions = self._parse_sacctmgr_output(self._run_command(f'sacctmgr --parsable2 show associations user={self.username}'))
        print(f"User Permissions: {self.user_permissions}")

    def _run_command(self, command):
        print(f"Running command: {command}")
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if stderr:
            print(f"Error executing command: {stderr.decode().strip()}")
        result = stdout.decode().strip()
        print(f"Command output: {result}")
        return result

    def _parse_scontrol_output(self, output):
        print("Parsing scontrol output")
        partitions_info = {}
        for line in output.split('\n'):
            parts = dict(re.findall(r'(\w+)=([^ ]+)', line))
            if 'PartitionName' in parts:
                partitions_info[parts['PartitionName']] = parts
        print(f"Parsed scontrol output: {partitions_info}")
        return partitions_info

    def _parse_sinfo_output(self, output):
        print("Parsing sinfo output")
        sinfo_data = {}
        for line in output.split('\n'):
            parts = line.split('|')
            if len(parts) == 3:
                partition_name = parts[0].strip('*')
                sinfo_data[partition_name] = {
                    'time_limit': parts[1],
                    'cpus': parts[2]
                }
        print(f"Parsed sinfo output: {sinfo_data}")
        return sinfo_data

    def _parse_sacctmgr_output(self, output):
        print("Parsing sacctmgr output")
        user_permissions = {}
        for line in output.split('\n'):
            parts = line.split('|')
            if len(parts) > 1 and parts[2] == self.username:
                user_permissions[parts[1]] = parts[-2].split(',')
        print(f"Parsed sacctmgr output: {user_permissions}")
        return user_permissions

    def find_best_partition_and_qos(self, runtime, cpus, gpus=0, gpu_type=None, constraints=None, account=None):
        if account is None:
            account = self._get_default_account()
        print(f"Using account: {account}")

        for partition, info in self.partitions_info.items():
            suitability, reason = self._is_partition_suitable(partition, runtime, cpus, gpus, gpu_type, constraints, account)
            print(f"Checking partition {partition}: Suitability: {suitability}, Reason: {reason}")
            if suitability:
                suitable_qos = self._get_best_qos(partition, account)
                print(f"Found suitable partition: {partition}, QoS: {suitable_qos}")
                return partition, suitable_qos

        print("No suitable partition found")
        return None, None

    def _is_partition_suitable(self, partition, runtime, cpus, gpus, gpu_type, constraints, account):
        partition_info = self.partitions_info.get(partition, {})
        
        if not self._is_user_or_group_allowed(partition_info, account):
            return False, "User or group not allowed"

        if not self._check_time_limit(partition, runtime):
            return False, "Time limit not met"

        if not self._check_cpu_availability(partition, cpus):
            return False, "CPU requirement not met"

        if gpus > 0 and not self._check_gpu_availability(partition_info, gpus):
            return False, "GPU requirement not met"

        return True, "Suitable"

    def _get_default_account(self):
        print("Fetching default account")
        command = f'sacctmgr --noheader --parsable2 show user {self.username} format=DefaultAccount'
        default_account = self._run_command(command).strip()
        print(f"Default account: {default_account}")
        return default_account    
    

    def _is_user_or_group_allowed(self, partition_info, account):
        allow_groups = partition_info.get('AllowGroups', '').split(',')
        deny_groups = partition_info.get('DenyGroups', '').split(',')
        allow_accounts = partition_info.get('AllowAccounts', '').split(',')
        deny_accounts = partition_info.get('DenyAccounts', '').split(',')

        if (self.username in allow_accounts) or (account in allow_accounts) or any(g in self.user_groups for g in allow_groups):
            if not any(g in self.user_groups for g in deny_groups) and self.username not in deny_accounts and account not in deny_accounts:
                return True, "User or group is allowed"
        return False, "User or group is not allowed"    
    

    def _check_time_limit(self, partition, job_time_limit):
        partition_time_limit = self.sinfo_data[partition]['time_limit']
        partition_minutes = self._time_limit_to_minutes(partition_time_limit)
        job_minutes = self._time_limit_to_minutes(job_time_limit)
        return job_minutes <= partition_minutes

    def _time_limit_to_minutes(self, time_limit):
        try:
            days, hours_minutes = time_limit.split('-', 1) if '-' in time_limit else (0, time_limit)
            hours, minutes = map(int, hours_minutes.split(':'))
            return int(days) * 1440 + hours * 60 + minutes
        except ValueError:
            return 0    
        
    def _check_cpu_availability(self, partition, job_cpus):
        total_cpus = int(self.sinfo_data[partition]['cpus'].split('/')[-1])
        return job_cpus <= total_cpus

    def _check_gpu_availability(self, partition_info, gpus):
        gres_str = partition_info.get('TRES', '')
        gpu_match = re.search(r'gres/gpu=(\d+)', gres_str)
        if gpu_match:
            available_gpus = int(gpu_match.group(1))
            return gpus <= available_gpus        
        return False        
    
    def _get_best_qos(self, partition, account):
        partition_qos = self.partitions_info.get(partition, {}).get('QoS', '').split(',')
        user_qos = self.user_permissions.get(account, [])
        common_qos = set(partition_qos).intersection(user_qos)
        return max(common_qos, default=None)


@click.command()
@click.option('--runtime', required=True, help='Requested runtime (e.g., "1-00:00:00" for 1 day).')
@click.option('--cpus', required=True, type=int, help='Number of CPUs requested.')
@click.option('--num_gpus', default=0, type=int, help='Number of GPUs requested.')
@click.option('--gpu_type', default=None, type=str, help='Type of GPU requested.')
@click.option('--constraints', default=None, type=str, help='Any additional constraints.')
@click.option('--account', default=None, type=str, help='Specific account to use.')

def cli(runtime, cpus, num_gpus, gpu_type, constraints, account):
    bb_class = BigBadClass()
    best_partition, best_qos = bb_class.find_best_partition_and_qos(runtime, cpus, num_gpus, gpu_type, constraints, account)
    print(f"Best Partition: {best_partition}, Best QoS: {best_qos}")

if __name__ == '__main__':
    cli()
