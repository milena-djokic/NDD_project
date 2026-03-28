#!/bin/bash
# Author: Milena Djokic

#SBATCH --job-name=ndd_project_af_interface_predictions
#SBATCH --partition=gpulong
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mincpus=64
#SBATCH --mem=512G
#SBATCH --time=7-00:00:00
#SBATCH --output=tmp1.log
#SBATCH --error=tmp2.log
#SBATCH --mail-user=m.djokic@imb-mainz.de
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT_50
#SBATCH --array=1-4425

#############################################################################################

# STAGE IN AND SCRATCH SETUP

# create the alphafold monitoring file:
export array_min=0
export job_name="ndd_project_af_interface_predictions"

job_id=${SLURM_JOB_ID}
array_job_id=${SLURM_ARRAY_JOB_ID}
task_id=${SLURM_ARRAY_TASK_ID}
pid=$$
echo "0" > "af_variable_${task_id}.txt"
# we wanna create some dissonance between the two scripts:
if (( $task_id == 1 )); then
	sleep 15
fi

if (( $task_id == 2 )); then
	sleep 2
fi

# Set up the input dir in the scratch space
export group_drive_dir="/fsimb/groups/imb-luckgr/imb-luckgr2/projects/AlphaFold/NDD_project_predictions/ndd_project_af_interface_predictions/"
cp -r "${group_drive_dir}input/" "$TMPDIR"
export INPUT_DIR=${TMPDIR}/input/
export INPUT_DIR_runs=${TMPDIR}/input/runs/

# Set up the output dir in the group drive:
export OUTPUT_DIR="${group_drive_dir}output/"
export OUTPUT_DIR_runs="${group_drive_dir}output/runs/"

if [ ! -d "${OUTPUT_DIR}" ]; then
	mkdir -p "$OUTPUT_DIR"
fi
if [ ! -d "${OUTPUT_DIR_runs}" ]; then
	mkdir -p "$OUTPUT_DIR_runs"
fi
#############################################################################################>

# Env setup - no need to load modules, just define a couple of variables
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_MEM_FRACTION=4.0

#############################################################################################>

CSV_FILE="alphafold_runs.tsv"
CSV_FILE_PATH="$INPUT_DIR$CSV_FILE"

prediction_names=()
run_ids=()
while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14 col15 col16 col17 col18 col19 col20 col21 col22 col23 col24 col25 col26
do
    # Skip the header line (optional)
    if [[ $col1 == "run_id" ]]; then
        continue
    fi

    # Print the variables
    prediction_names+=($col6)
    run_ids+=("$col1")
done < "$CSV_FILE_PATH"

###########################################

export running_dictionary_file_path="running_dictionary.txt"
export finalized_dictionary_file_path="finalized_dictionary.txt"

touch "$running_dictionary_file_path"
touch "$finalized_dictionary_file_path"
touch "first_time_running_MSA.txt"
touch "recycling_MSA.txt"
touch "repeating_MSA.txt"
touch "malfunctioning_MSA_tracking.txt"
touch "all_frag.txt"

if [ ! -s "$running_dictionary_file_path" ] && [ $array_min -eq $task_id ]; then
	# get a list of fragments and let's store this information of which fragment already ran and which did not yet
	temp_file=$(mktemp) # this creates a unique file
	for item in "${prediction_names[@]}"; do
		fragment="${item#*_}"
		echo "$fragment" | tr '.' '\n'
	done | sort | uniq > "$temp_file"
	IFS=$'\n' read -r -d '' -a fragments < "$temp_file" # read the fragments into an array
	rm "$temp_file"

	for fragment in "${fragments[@]}"; do
		echo "$fragment 0" >> $running_dictionary_file_path
	done
	for fragment in "${fragments[@]}"; do
		echo "$fragment 0" >> $finalized_dictionary_file_path
	done
fi

# Run the monitoring script in the background!
source ./monitoring_script.sh ${prediction_names[$SLURM_ARRAY_TASK_ID]} $job_id $array_job_id $task_id $pid &
sleep 1

source ./single_job_script.sh ${prediction_names[$SLURM_ARRAY_TASK_ID]} $job_id $array_job_id $task_id $pid \
1>> stdout.log 2>> stderr.log

# STAGE OUT:

# and just to make sure that all MSAs are transferred into the central directory
unique_run_ids=($(echo "${run_ids[@]}" | tr ' ' '\n' | sort -n | uniq))
for run_id in "${unique_run_ids[@]}"; do
    fasta_path="${INPUT_DIR_runs}${run_id}/" # input dir from the scratch directory
    run_path="${OUTPUT_DIR_runs}run${run_id}/" # output dir on the group drive
    store_path="${group_drive_dir}output/runs/run${run_id}/" # actually now same as the run path
    msas_path="${group_drive_dir}output/computed_msas/" # central directory for recycling msas
    # and these dirs are not run/prediction specific but rather for the whole protein pair
    prediction_name="None"
    python3 /fsimb/groups/imb-luckgr/imb-luckgr2/projects/AlphaFold/NDD_project_predictions/recycle_precomputed_msas.py \
                                      --prediction_name $prediction_name \
                                      --fasta_path $fasta_path \
                                      --path_to_run $run_path \
                                      --storage_path $store_path \
                                      --precomputed_msas_path $msas_path
done

# remove temporary files
rm tmp1.log
rm tmp2.log







