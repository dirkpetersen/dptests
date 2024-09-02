#! /usr/bin/env python3

import subprocess
import os
import re
import grp

class BigBadClass:
    def __init__(self):
        self.username = os.getenv("USER")
        self.user_groups = [g.gr_name for g in grp.getgrall() if self.username in g.gr_mem]
        self.partitions_info = self._parse_scontrol_output(self._run_command('scontrol --oneline show partition'))
        self.sinfo_data = self._parse_sinfo_output(self._run_command('sinfo -o "%P|%l|%C"'))
        self.user_permissions = self._parse_sacctmgr_output(self._run_command(f'sacctmgr -p show associations user={self.username}'))

    def _run_command(self, command):
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if stderr:
            raise Exception(f"Error executing command: {stderr.decode().strip()}")
        return stdout.decode().strip()

    def _parse_scontrol_output(self, output):
        partitions_info = {}
        for line in output.split('\n'):
            parts = dict(re.findall(r'(\w+)=([^ ]+)', line))
            if 'PartitionName' in parts:
                partitions_info[parts['PartitionName']] = parts
        return partitions_info

    def _parse_sinfo_output(self, output):
        sinfo_data = {}
        for line in output.split('\n'):
            parts = line.split('|')
            if len(parts) == 3:
                partition_name = parts[0].strip('*')
                sinfo_data[partition_name] = {
                    'time_limit': parts[1],
                    'cpus': parts[2]
                }
        return sinfo_data

    def _parse_sacctmgr_output(self, output):
        user_permissions = {}
        for line in output.split('\n'):
            parts = line.split('|')
            if len(parts) > 1 and parts[2] == self.username:
                user_permissions[parts[1]] = parts[-2].split(',')
        return user_permissions

    def find_best_partition_and_qos(self, runtime, cpus, num_gpus=0, gpu_type=None, constraints=None, account=None):
        if account is None:
            account = self._get_default_account()

        for partition, info in self.partitions_info.items():
            suitability, reason = self._is_partition_suitable(partition, runtime, cpus, num_gpus, gpu_type, constraints, account)
            if suitability:
                suitable_qos = self._get_best_qos(partition, account)
                if suitable_qos:
                    return partition, suitable_qos
            else:
                print(f"Partition {partition} is unsuitable: {reason}")
        return None, None

    def _is_partition_suitable(self, partition, runtime, cpus, num_gpus, gpu_type, constraints, account):
        partition_info = self.partitions_info.get(partition, {})
        
        if not self._is_user_or_group_allowed(partition_info, account):
            return False, "User or group not allowed"

        if not self._check_time_limit(partition, runtime):
            return False, "Time limit not met"

        if not self._check_cpu_availability(partition, cpus):
            return False, "CPU requirement not met"

        if num_gpus > 0 and not self._check_gpu_availability(partition_info, num_gpus):
            return False, "GPU requirement not met"

        return True, "Suitable"

    def _get_best_qos(self, partition, account):
        partition_qos = self.partitions_info.get(partition, {}).get('QoS', '').split(',')
        user_qos = self.user_permissions.get(account, [])
        common_qos = set(partition_qos).intersection(user_qos)
        return max(common_qos, default=None)

    def _is_user_or_group_allowed(self, partition_info, account):
        allow_groups = partition_info.get('AllowGroups', '').split(',')
        deny_groups = partition_info.get('DenyGroups', '').split(',')
        allow_accounts = partition_info.get('AllowAccounts', '').split(',')
        deny_accounts = partition_info.get('DenyAccounts', '').split(',')

        if any(group in self.user_groups for group in allow_groups) or self.username in allow_accounts or account in allow_accounts:
            return True
        if any(group in self.user_groups for group in deny_groups) or self.username in deny_accounts or account in deny_accounts:
            return False
        return True

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

    def _check_gpu_availability(self, partition_info, num_gpus):
        gres_str = partition_info.get('TRES', '')
        gpu_match = re.search(r'gres/gpu=(\d+)', gres_str)
        if gpu_match:
            available_gpus = int(gpu_match.group(1))
            return num_gpus <= available_gpus
        return False

    def _get_default_account(self):
        command = 'sacctmgr --noheader --parsable2 show user {} format=DefaultAccount'.format(self.username)
        return self._run_command(command).strip()

# Example usage
bb_class = BigBadClass()
best_partition, best_qos = bb_class.find_best_partition_and_qos(runtime='1-00:00:00', cpus=16, num_gpus=2, gpu_type='Tesla', constraints='high-mem')
print(f"Best Partition: {best_partition}, Best QoS: {best_qos}")


