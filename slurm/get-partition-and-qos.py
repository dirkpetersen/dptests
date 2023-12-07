
import subprocess
import grp
import os

def get_shell_output(command):
    """Execute a shell command and return its output."""
    result = subprocess.run(command, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return result.stdout.strip()

def parse_tabular_data(data_str, separator="|"):
    """Parse data (e.g. acctmgr) presented in a tabular format into a list of dictionaries."""
    lines = data_str.strip().splitlines()
    headers = lines[0].split(separator)
    data = []
    for line in lines[1:]:
        values = line.split(separator)
        data.append(dict(zip(headers, values)))
    return data

def parse_partition_data(data_str):
    """Parse data presented in a tabular format into a list of dictionaries."""
    lines = data_str.strip().split('\n')
    # Parse each line into a dictionary
    partitions = []
    for line in lines:
        parts = line.split()
        partition_dict = {}
        for part in parts:
            key, value = part.split("=", 1)
            partition_dict[key] = value
        partitions.append(partition_dict)
    return partitions

def get_user_groups():
    """Get the groups the current Unix user is a member of."""
    groups = [grp.getgrgid(gid).gr_name for gid in os.getgroups()]
    return groups

def convert_to_hours(time_str):
    """Convert a time string to hours."""
    days, hours, minutes = 0, 0, 0
    if '-' in time_str:
        days, time_str = time_str.split('-')
        days = int(days) * 24
    if ':' in time_str:
        hours, minutes = map(int, time_str.split(':'))
    return days + hours + minutes / 60

def find_best_qos(allowed_qos, qos_info, required_runtime_hours):
    """Find the best QoS based on runtime requirements."""
    best_qos = None
    shortest_runtime = float('inf')

    for qos in allowed_qos:
        max_wall = convert_to_hours(qos_info.get(qos, {}).get('MaxWall', '0'))
        if max_wall >= required_runtime_hours and max_wall < shortest_runtime:
            best_qos = qos
            shortest_runtime = max_wall

    return best_qos

def get_allowed_partitions_and_qos(partitions, associations, user_groups, account=None):
    """Get a dictionary with keys = partitions and values = QOSs the user is allowed to use."""
    bacc=os.environ.get('SBATCH_ACCOUNT', '')
    account = bacc if bacc else account
    sacc=os.environ.get('SLURM_ACCOUNT', '')
    account = sacc if sacc else account
    allowed_partitions = {}
    for partition in partitions:
        pname=partition['PartitionName']
        add_partition = False
        if partition.get('State', '') != 'UP':
            continue
        if any(group in user_groups for group in partition.get('DenyGroups', '').split(',')):
            continue
        if account in partition.get('DenyAccounts', '').split(','):            
            continue
        allowedaccounts=partition.get('AllowAccounts', '').split(',')
        if allowedaccounts != ['']:
            if account in allowedaccounts or 'ALL' in allowedaccounts:
                add_partition = True
        elif any(group in user_groups for group in partition.get('AllowGroups', '').split(',')):
            add_partition = True
        elif partition.get('AllowGroups', '') == 'ALL':
            add_partition = True
        if add_partition:
            p_deniedqos = partition.get('DenyQos', '').split(',')            
            p_allowedqos = partition.get('AllowQos', '').split(',')            
            account_qos = associations.get(account, [])
            if p_deniedqos != ['']:
                allowed_qos = [q for q in account_qos if q not in p_deniedqos]
                #print(f"p_deniedqos: allowed_qos in {pname}:", allowed_qos) 
            elif p_allowedqos == ['ALL']:
                allowed_qos = account_qos
                #print(f"p_allowedqos = ALL in {pname}:", allowed_qos)                 
            elif p_allowedqos != ['']:
                allowed_qos = [q for q in account_qos if q in p_allowedqos]
                #print(f"p_allowedqos: allowed_qos in {pname}:", allowed_qos) 
            else:
                allowed_qos = []
                #print(f"p_allowedqos = [] in {pname}:", allowed_qos)     
            allowed_partitions[pname] = allowed_qos
    return allowed_partitions

def find_best_partition_and_qos(num_cores, runtime_hours, partitions, associations, user_groups, qos_info):
    """Find the best partition and QoS based on given criteria."""
    best_partition = None
    best_qos = None

    for partition in partitions:
        # Check user group and account access
        if (any(group in user_groups for group in partition.get('AllowGroups', '').split(',')) and
            not any(group in user_groups for group in partition.get('DenyGroups', '').split(',')) and
            any(account in associations for account in partition.get('AllowAccounts', '').split(',')) and
            not any(account in associations for account in partition.get('DenyAccounts', '').split(',')) and
            int(partition.get('TotalCPUs', 0)) >= num_cores and
            float(partition.get('MaxTime', '0').replace('-', '')) >= runtime_hours):

            best_partition = partition['PartitionName']
            allowed_qos = partition.get('AllowQos', '').split(',')
            best_qos = find_best_qos(allowed_qos, qos_info, runtime_hours)
            if best_qos:
                break

    return best_partition, best_qos

# Retrieve partition information
partition_str = get_shell_output("scontrol show partition --oneliner")
partitions = parse_partition_data(partition_str)

#print('Partitions:', partitions)

# Retrieve associations information for the current user
current_user = os.getlogin()
print('Current User:', current_user)
associations_str = get_shell_output(f"sacctmgr show associations where user={current_user} format=Account,QOS --parsable2")
associations = {item['Account']: item['QOS'].split(",") for item in parse_tabular_data(associations_str) if 'Account' in item}
print('Associations:', associations)

default_account = get_shell_output(f'sacctmgr --noheader --parsable2 show user {current_user} format=DefaultAccount')

# Retrieve QoS information
qos_str = get_shell_output("sacctmgr show qos --parsable2")
qos_info = {item['Name']: item for item in parse_tabular_data(qos_str)}
#print('QOS Info:', qos_info)

# Get the Unix groups for the current user
user_groups = get_user_groups()
print('User Groups:', user_groups)

allowed_partitons=get_allowed_partitions_and_qos(partitions, associations, user_groups, default_account)
#print('Allowed Partitions:', allowed_partitons)

print('***** Allowed Partitions:')
for k, v in allowed_partitons.items():
    print(k,v)

# Find best partition and QoS
num_cores = 10
runtime_hours = 20
best_partition, best_qos = find_best_partition_and_qos(num_cores, runtime_hours, partitions, associations, user_groups, qos_info)
print(f"Best Partition: {best_partition}")
print(f"Best QoS: {best_qos}")

