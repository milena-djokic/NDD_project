#!/bin/bash
# Author: Milena Djokic
#############################################################################################
# Load the csv file and all the variables in:
CSV_FILE_PATH="${job_name}/input/alphafold_runs.tsv"

run_id=()
run_task_id=()
task_id=()
fasta_file=()
output_dir=()
log_dir=()
usage_dir=()
monitoring_dir=()

prediction_name=$1
job_id_monitoring=$2
array_job_id_monitoring=$3
task_id_monitoring=$4
pid=$5
array_job_task_id_monitoring="${array_job_id_monitoring}_${task_id_monitoring}"

# Read the data we need from the csv file
while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14 col15 col16 col17 col18 col19 col20 col21 col22 col23 col24 col25 col26
do
    # Skip the header line (optional)
    if [[ $col6 == $prediction_name ]]; then
        # Print the variables
        run_id=$col1
        run_task_id=$col2
        task_id=$col3
        fasta_file=$col22
        output_dir=$col23
        log_dir=$col24
        usage_dir=$col25
        monitoring_dir=$col26
    else
        continue
    fi
done < "$CSV_FILE_PATH"

#############################################################################################
# let's first record the allocated resources:

# number of CPUs
# number of cores per CPU - its actually all cores so we don't need to track this for now
# number of GPUs
# number of cores per GPU - its also all cores, so we don't need to track that either
# memory

 # we need to specify the exact path later and then we can also remove the job task id
alloc_resources=$(mktemp)
scontrol show job $array_job_task_id_monitoring >> "$alloc_resources"

job_state=$(grep -oP 'JobState=\K\w+' $alloc_resources)
if [[ $job_state == "RUNNING" ]]; then #to make sure we got the correct job and not some that are no longer running
	# fetch the TRACKABLE RESOURCES (TRES):
	alloc_cpus=$(grep 'TRES=' $alloc_resources | grep -o "cpu=[0-9]*" | cut -d= -f2)
	alloc_gpus=$(grep 'TRES=' $alloc_resources | grep -o "gres/gpu=[0-9]*" | cut -d= -f2)
	alloc_mem=$(grep 'TRES=' $alloc_resources | grep -o "mem=[0-9]*" | cut -d= -f2)
	alloc_mem="${alloc_mem} G"
fi

rm "$alloc_resources"
# and now the dynamic ones I should probably wrap around in a function, I am not sure, we will see about this later

#############################################################################################

# grep current usage of trackable resources:

track_resources() {

	header=$(echo "$(sstat -ap -j "$job_id_monitoring")" | head -n 1)
	batch_line=$(echo "$(sstat -ap -j "$job_id_monitoring")" | grep "$job_id_monitoring\.batch")

	# Split the header and batch_line by pipes (|) into arrays
	IFS='|' read -r -a header_fields <<< "$header"
	IFS='|' read -r -a batch_fields <<< "$batch_line"

	# Function to get the position of a field in the header
	get_field_position() {
	  local field_name="$1"
	  for i in "${!header_fields[@]}"; do
	    if [[ "${header_fields[$i]}" == "$field_name" ]]; then
	      echo "$i"
	      return
	    fi
	  done
	}

	# Get the dynamic positions of the TRES fields
	tres_in_ave_pos=$(get_field_position "TRESUsageInAve")
	tres_in_max_pos=$(get_field_position "TRESUsageInMax")
	tres_in_min_pos=$(get_field_position "TRESUsageInMin")
	tres_in_tot_pos=$(get_field_position "TRESUsageInTot")

	# Extract the corresponding TRES values from the batch line using the positions
	tres_usage_in_ave=${batch_fields[$tres_in_ave_pos]}
	tres_usage_in_max=${batch_fields[$tres_in_max_pos]}
	tres_usage_in_min=${batch_fields[$tres_in_min_pos]}
	tres_usage_in_tot=${batch_fields[$tres_in_tot_pos]}

	cpu_in_ave=$(echo "$tres_usage_in_ave" | grep -o "cpu=[^,]*" | cut -d= -f2)
	mem_in_ave=$(echo "$tres_usage_in_ave" | grep -o ",mem=[^,]*" | cut -d= -f2)

	cpu_in_max=$(echo "$tres_usage_in_max" | grep -o "cpu=[^,]*" | cut -d= -f2)
	mem_in_max=$(echo "$tres_usage_in_max" | grep -o ",mem=[^,]*" | cut -d= -f2)

	cpu_in_min=$(echo "$tres_usage_in_min" | grep -o "cpu=[^,]*" | cut -d= -f2)
	mem_in_min=$(echo "$tres_usage_in_min" | grep -o ",mem=[^,]*" | cut -d= -f2)

	cpu_in_tot=$(echo "$tres_usage_in_tot" | grep -o "cpu=[^,]*" | cut -d= -f2)
	mem_in_tot=$(echo "$tres_usage_in_tot" | grep -o ",mem=[^,]*" | cut -d= -f2)

	# change the format of the mem variables!
	mem_in_ave="${mem_in_ave//K/ K}"
	mem_in_max="${mem_in_max//K/ K}"
	mem_in_min="${mem_in_min//K/ K}"
	mem_in_tot="${mem_in_tot//K/ K}"

	# calculate the CPU usage:
	# Function to convert HH:MM:SS to seconds
	convert_to_seconds() {
	  local time=$1
	  local hours=$(echo "$time" | cut -d':' -f1)
	  local minutes=$(echo "$time" | cut -d':' -f2)
	  local seconds=$(echo "$time" | cut -d':' -f3)
	  # Remove leading zeros
	  # remove the leading 0
	  hours="${hours#0}"
	  minutes="${minutes#0}"
	  seconds="${seconds#0}"
	  # Handle the case where the entire number was '0'
	  hours="${hours:-0}"
	  minutes="${minutes:-0}"
	  seconds="${seconds:-0}"
	  echo $((hours * 3600 + minutes * 60 + seconds))
	}

	# Ensure function is called properly and variables are correctly assigned
	cpu_in_ave_seconds=$(convert_to_seconds "$cpu_in_ave")
	cpu_in_max_seconds=$(convert_to_seconds "$cpu_in_max")
	cpu_in_min_seconds=$(convert_to_seconds "$cpu_in_min")
	cpu_in_tot_seconds=$(convert_to_seconds "$cpu_in_tot")

	runtime=$(echo "$(scontrol show job "$array_job_task_id_monitoring")" | grep -oP 'RunTime=\K[^ ]+')
	runtime_seconds=$(convert_to_seconds "$runtime")
	# Check if the value is 0
	if [ "$runtime_seconds" -eq 0 ]; then
		cpu_in_ave_usage="None"
		cpu_in_max_usage="None"
		cpu_in_min_usage="None"
		cpu_in_tot_usage="None"
	else
		cpu_in_ave_usage=$(echo "scale=2; ($cpu_in_ave_seconds / $runtime_seconds) * 100" | bc)
		cpu_in_max_usage=$(echo "scale=2; ($cpu_in_max_seconds / $runtime_seconds) * 100" | bc)
		cpu_in_min_usage=$(echo "scale=2; ($cpu_in_min_seconds / $runtime_seconds) * 100" | bc)
		cpu_in_tot_usage=$(echo "scale=2; ($cpu_in_tot_seconds / $runtime_seconds) * 100" | bc)

		# reformat those
		cpu_in_ave_usage=$(printf "%.0f" "$cpu_in_ave_usage")
		cpu_in_ave_usage="$cpu_in_ave_usage %"
		cpu_in_max_usage=$(printf "%.0f" "$cpu_in_max_usage")
		cpu_in_max_usage="$cpu_in_max_usage %"
		cpu_in_min_usage=$(printf "%.0f" "$cpu_in_min_usage")
		cpu_in_min_usage="$cpu_in_min_usage %"
		cpu_in_tot_usage=$(printf "%.0f" "$cpu_in_tot_usage")
		cpu_in_tot_usage="$cpu_in_tot_usage %"
	fi



	# to see if alphafold is currently running or not
	if [ -f "af_variable_${task_id_monitoring}.txt" ]; then
		alphafold_running=$(cat "af_variable_${task_id_monitoring}.txt")
	fi

	output=$(srun --jobid="$job_id_monitoring" nvidia-smi --query-gpu=timestamp,gpu_uuid,index,utilization.gpu,pstate,compute_mode,temperature.gpu,power.draw,power.draw.average,enforced.power.limit,memory.total,memory.reserved,memory.used,memory.free,utilization.memory --format=csv,noheader)
	# Set IFS to comma
	IFS=',' read -r timestamp gpu_uuid index utilization_gpu pstate compute_mode temperature_gpu power_draw power_draw_average enforced_power_limit memory_total memory_reserved memory_used memory_free utilization_memory <<< "$output"

	# Trim whitespace
	timestamp=$(echo $timestamp | xargs)
	gpu_uuid=$(echo $gpu_uuid | xargs)
	index=$(echo $index | xargs)
	utilization_gpu=$(echo $utilization_gpu | xargs)
	pstate=$(echo $pstate | xargs)
	compute_mode=$(echo $compute_mode | xargs)
	temperature_gpu=$(echo $temperature_gpu | xargs)
	power_draw=$(echo $power_draw | xargs)
	power_draw_average=$(echo $power_draw_average | xargs)
	enforced_power_limit=$(echo $enforced_power_limit | xargs)
	memory_total=$(echo $memory_total | xargs)
	memory_reserved=$(echo $memory_reserved | xargs)
	memory_used=$(echo $memory_used | xargs)
	memory_free=$(echo $memory_free | xargs)
	utilization_memory=$(echo $utilization_memory | xargs)

	#reformat some of these variables
	timestamp=$(echo "$timestamp" | sed 's/\//-/g; s/\.[0-9]\{1,\}//')
	memory_total=$(echo "$memory_total" | sed 's/ MiB/ M/')
	memory_reserved=$(echo "$memory_reserved" | sed 's/ MiB/ M/')
	memory_used=$(echo "$memory_used" | sed 's/ MiB/ M/')
	memory_free=$(echo "$memory_free" | sed 's/ MiB/ M/')

	# change column names and order
	alphafold_running="$alphafold_running"

	cpu_n_allocated="$alloc_cpus"
	cpu_seconds_ave="$cpu_in_ave"
	cpu_seconds_max="$cpu_in_max"
	cpu_seconds_min="$cpu_in_min"
	cpu_seconds_tot="$cpu_in_tot"
	cpu_usage_ave="$cpu_in_ave_usage"
	cpu_usage_max="$cpu_in_max_usage"
	cpu_usage_min="$cpu_in_min_usage"
	cpu_usage_tot="$cpu_in_tot_usage"
	cpu_ps_usage=$(ps -p $pid -o %cpu | tail -n 1) # $$ this gives us PID, but this will change so idk
	cpu_ps_usage="${cpu_ps_usage} %"
	cpu_ps_n_cpus=$(ps -o nlwp= -p $pid | grep -o '[0-9]*')

	gpu_n_allocated="$alloc_gpus"
	gpu_id="$gpu_uuid"
	gpu_index="$index"
	gpu_percent_utilization="$utilization_gpu"
	gpu_performance_state="$pstate"
	gpu_compute_mode="$compute_mode"
	gpu_temperature="$temperature_gpu"
	gpu_power_draw="$power_draw"
	gpu_power_draw_average="$power_draw_average"
	gpu_power_limit="$enforced_power_limit"
	gpu_total_memory="$memory_total"
	gpu_resrved_memory="$memory_reserved"
	gpu_used_memory="$memory_used"
	gpu_free_memory="$memory_free"
	gpu_memory_percent_utilization="$utilization_memory"

	memory_allocated="$alloc_mem"
	memory_ave="$mem_in_ave"
	memory_max="$mem_in_max"
	memory_min="$mem_in_min"
	memory_tot="$mem_in_tot"
	
	# we do this in the end so that it's as accurate as possible
	currtime="$(date '+%Y-%m-%d %H:%M:%S')"
	runtime=$(echo "$(scontrol show job "$array_job_task_id_monitoring")" | grep -oP 'RunTime=\K[^ ]+')

	timestamp="$currtime" # this one and our manual ones are the same so I will omit one of them
	runtime="$runtime"

	# and now push all of that to a tsv file
	tsv_file="${OUTPUT_DIR_runs}run${run_id}/run${prediction_name}/log_dir/usage_dir/resource_usage_monitoring.tsv"
	touch $tsv_file

	if [ ! -s "$tsv_file" ]; then # if the file is empty, then add a header with column names
	    echo -e "timestamp\truntime\talphafold_running\tcpu_n_allocated\tcpu_seconds_ave\tcpu_seconds_max\tcpu_seconds_min\tcpu_seconds_tot\tcpu_usage_ave\tcpu_usage_max\tcpu_usage_min\tcpu_usage_tot\tcpu_ps_usage\tcpu_ps_n_cpus\tgpu_n_allocated\tgpu_id\tgpu_index\tgpu_percent_utilization\tgpu_performance_state\tgpu_compute_mode\tgpu_temperature\tgpu_power_draw\tgpu_power_draw_average\tgpu_power_limit\tgpu_total_memory\tgpu_resrved_memory\tgpu_used_memory\tgpu_free_memory\tgpu_memory_percent_utilization\tmemory_allocated\tmemory_ave\tmemory_max\tmemory_min\tmemory_tot" > "$tsv_file"
	fi

	# Append variables as a single row in a TSV format
	echo -e "$timestamp\t$runtime\t$alphafold_running\t$cpu_n_allocated\t$cpu_seconds_ave\t$cpu_seconds_max\t$cpu_seconds_min\t$cpu_seconds_tot\t$cpu_usage_ave\t$cpu_usage_max\t$cpu_usage_min\t$cpu_usage_tot\t$cpu_ps_usage\t$cpu_ps_n_cpus\t$gpu_n_allocated\t$gpu_id\t$gpu_index\t$gpu_percent_utilization\t$gpu_performance_state\t$gpu_compute_mode\t$gpu_temperature\t$gpu_power_draw\t$gpu_power_draw_average\t$gpu_power_limit\t$gpu_total_memory\t$gpu_resrved_memory\t$gpu_used_memory\t$gpu_free_memory\t$gpu_memory_percent_utilization\t$memory_allocated\t$memory_ave\t$memory_max\t$memory_min\t$memory_tot" >> $tsv_file

}

while true; do
	track_resources
	sleep 1
done








