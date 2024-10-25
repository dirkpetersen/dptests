#! /bin/bash

# Constants for column headers
MEMORY_UTIL_COL="MEMORY_UTIL"
PROCESS_NAME_COL="PROCESS_NAME"
USER_COL="USER"
DATE_TIME_COL="DATE_TIME"

# Display header with extra width for User field
printf "%-20s %-6s %-20s %-8s %-15s %-15s %-15s\n" "${DATE_TIME_COL}" "GPU ID" \
       "${USER_COL}" "PID" "${MEMORY_UTIL_COL}" "NVLINK_UTIL" "${PROCESS_NAME_COL}"

# Function to check GPU usage, memory utilization, and NVLink status for non-root users
monitor_gpu_vram_per_user() {
  for gpu_id in $(nvidia-smi --query-gpu=index --format=csv,noheader); do

    # Determine NVLink status
    nvlink_status=$(nvidia-smi nvlink --status --id=${gpu_id} 2>&1)
    if echo "${nvlink_status}" | grep -q "inActive"; then
      nvlink_usage="inactive"
    elif echo "${nvlink_status}" | grep -q "Unable to retrieve NVLink information"; then
      nvlink_usage="N/A"
    else
      nvlink_usage=$(echo "${nvlink_status}" | grep 'Link' | awk '{print $3 " " $4 " " $5}')
    fi

    # Capture all processes using the GPU
    nvidia-smi --query-compute-apps=pid,used_memory,process_name \
      --id=${gpu_id} --format=csv,noheader,nounits \
      | while IFS=, read -r pid vram_usage process_name; do
        # Fetch user of the process
        user=$(ps -o user= -p "${pid}" 2>/dev/null | awk 'NR==1{print $1}')

        # Skip root processes
        if [[ ${user} != "root" ]]; then
          # Get current date and time
          datetime=$(date "+%Y-%m-%d %H:%M:%S")

          # Output results for this process with date and time, wider User field, and NVLink status
          printf "%-20s %-6s %-20s %-8s %-15s %-15s %-15s\n" \
                 "${datetime}" "${gpu_id}" "${user}" "${pid}" "${vram_usage} MB" \
                 "${nvlink_usage}" "${process_name}"
        fi
      done
  done
}

# Run the function
monitor_gpu_vram_per_user

