import subprocess
import os

# This creates multiple Virtual Box hosts to build a PetaSAN cluster 
# before running this script, execute this command 4 times to create 4 hostonly NICS:
#    VBoxManage hostonlyif create
# Configure the NICs in the VirtualBox with the IP space you want to use for management, iscsi, cifs, s3 etc 
# e.g. 10.0.1.0/24, 10.0.2.0/24, 10.0.3.0/24, 10.0.4.0/24  

# Configuration
vboxmanage = "C:/Program Files/Oracle/VirtualBox/VBoxManage.exe"
iso_path = "D:/ISO/petasan-3.2.0.iso"
vmdir = "D:/Users/dipei/VirtualBox"
nicname = "VirtualBox Host-Only Ethernet Adapter"
#nicbridge = "Intel(R) Wi-Fi 6E AX211 160MHz"  # no longer needed
cpus = 2
memory = 4096  # in MB
disk_size = 65536  # 64GB or more, in MB
machine_count = 3
network_card_count = 4
disk_count = 8

if not os.path.exists(vboxmanage):
    vboxmanage = "VBoxManage"

for j in range(machine_count):
    # VM configuration
    vm_name = f"PetaSAN_{j+1}"

    # Create the VM
    subprocess.run([vboxmanage, "createvm", "--name", vm_name, "--ostype", "SUSE_LE_64", "--register"], check=True)

    os.chdir(os.path.join(vmdir,vm_name))

    # Set the VM's CPU and memory
    subprocess.run([vboxmanage, "modifyvm", vm_name, "--cpus", str(cpus), "--memory", str(memory), "--graphicscontroller=vmsvga"],check=True)

    # Attach the SUSE ISO to the VM as a DVD drive
    subprocess.run([vboxmanage, "storagectl", vm_name, "--name", "IDE Controller", "--add", "ide"], check=True)
    subprocess.run([vboxmanage, "storageattach", vm_name, "--storagectl", "IDE Controller", "--port", "0", "--device", "0", "--type", "dvddrive", "--medium", iso_path], check=True)

    # Create and attach the virtual hard disks
    subprocess.run([vboxmanage, "storagectl", vm_name, "--name", "SAS Controller", "--add", "sas", "--controller", "LSILogicSAS", "--portcount", str(disk_count+1)], check=True)
    for i in range(disk_count+1):
        disk_path = f"{vm_name}_disk_{i}.vdi"
        subprocess.run([vboxmanage, "createhd", "--filename", disk_path, "--size", str(disk_size)], check=True)
        subprocess.run([vboxmanage, "storageattach", vm_name, "--storagectl", "SAS Controller", "--port", str(i), "--device", "0", "--type", "hdd", "--medium", disk_path], check=True)

    # Add network cards
    for i in range(network_card_count):
        mynic=nicname
        if i > 0:
            mynic = f"{nicname} #{i+1}"
        subprocess.run([vboxmanage, "modifyvm", vm_name, "--nic" + str(i + 1), "hostonly", "--hostonlyadapter" + str(i + 1), mynic], check=True) 
        #subprocess.run([vboxmanage, "modifyvm", vm_name, "--nic" + str(i + 1), "bridged", "--bridge-adapter" + str(i + 1), nicbridge], check=True)
