#!/bin/bash
# Author: Milena Djokic

# Set up variables
MODEL_PRESET="multimer"
# Relaxation on GPU
GPU_RELAX="true"

# initiate the input variables
prediction_name=$1
job_id_here=$2
array_job_id_here=$3
task_id_here=$4
pid=$5
job_task_id_here="${array_job_id_here}_${task_id_here}"

#####################################################################################

# Check if we've already done predictions for one of these two fragments - load into variables
fragments="${prediction_name#[0-9]_}"
fragment1="${fragments%%.*}"
fragment2="${fragments#*.}"

frag1_started=$(awk -v frag="$fragment1" '$1 == frag {print $2}' $running_dictionary_file_path)
frag1_finalized=$(awk -v frag="$fragment1" '$1 == frag {print $2}' $finalized_dictionary_file_path)
frag2_started=$(awk -v frag="$fragment2" '$1 == frag {print $2}' $running_dictionary_file_path)
frag2_finalized=$(awk -v frag="$fragment2" '$1 == frag {print $2}' $finalized_dictionary_file_path)
# the values are 0 for all 4 here, but we have to create some conditions
fragA="${frag1_started}_${frag1_finalized}"
fragB="${frag2_started}_${frag2_finalized}"

# Change the value to 1 as soon as possible, but we have to lock the file to make sure that it is accessed by one job at the time!
(
  flock -x 200  # Exclusive lock on file descriptor 200

  new_val=1
  temp_file=$(mktemp) # this creates a unique file
  awk -v frag1="$fragment1" -v frag2="$fragment2" -v new_val="$new_val" '
  {
    if ($1 == frag1 || $1 == frag2) {
      # Replace the value for either fragment1 or fragment2
      print $1, new_val
    } else {
      # Print unchanged lines
      print
    }
  }' "$running_dictionary_file_path" > "$temp_file" && mv "$temp_file" "$running_dictionary_file_path"

) 200>"$running_dictionary_file_path.lock"


####################################################################################

# We are running this code in the monitoring script:
# Load the csv file and all the variables in:
CSV_FILE_PATH="$TMPDIR/input/alphafold_runs.tsv"

run_id=()
run_task_id=()
task_id=()
fasta_file=()
output_dir=()
log_dir=()
usage_dir=()
monitoring_dir=()

# Read the data we need from the csv file
while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14 col15 col16 col17 col18 col19 col20 col21 col22 col23 col24 col25 col26
do
    # Skip the header line (optional)
    if [[ $col6 == $prediction_name ]]; then
        # Print the variables
        run_id=$col1 #
        run_task_id=$col2
        task_id=$col3
        fasta_file=$col22
        output_dir=$col23 #
        log_dir=$col24 #
        usage_dir=$col25 #
        monitoring_dir=$col26 #
    else
        continue
    fi
done < "$CSV_FILE_PATH"

# Print out the current process name, as well as which GPU is being used by the script
task_id_here_1=$((task_id_here + 1))
echo "${run_id} ${run_task_id} ${task_id_here_1} ${prediction_name} ${SLURM_JOB_GPUS}"

# Create the run and the prediction directory
if [ ! -d "$OUTPUT_DIR_runs$run_id" ]; then
    mkdir -p "${OUTPUT_DIR_runs}run${run_id}/"
    mkdir -p "${OUTPUT_DIR_runs}${output_dir}log_dir/usage_dir/"
    mkdir -p "${OUTPUT_DIR_runs}${output_dir}log_dir/monitoring_dir/"
else # otherwise just create the prediction directory, unless it already exists!
    if [ ! -d "${OUTPUT_DIR_runs}${output_dir}log_dir/usage_dir/" ]; then
        mkdir -p "${OUTPUT_DIR_runs}${output_dir}log_dir/usage_dir/"
        mkdir -p "${OUTPUT_DIR_runs}${output_dir}log_dir/monitoring_dir/"
    fi
fi

run_dir="${OUTPUT_DIR_runs}${output_dir}"
fasta_file_path="${INPUT_DIR_runs}${fasta_file}"

####################################################################################

# Now we run the AF predictions, depending on whether we ran this before or not:

if [[ $fragA == "1_0" || $fragA == "1_1" || $fragB == "1_0" || $fragB == "1_1" ]]; then
	# let's try to recycle the MSA either way, maybe 1_0 has finished by now!
	fasta_path="${INPUT_DIR_runs}${run_id}/" # input dir from the scratch directory
	run_path="${OUTPUT_DIR_runs}run${run_id}/" # output dir on the group drive
	store_path="${group_drive_dir}output/runs/run${run_id}/" # actually now same as the run path
	msas_path="${group_drive_dir}output/computed_msas/" # central directory for recycling msas
	# and these dirs are not run/prediction specific but rather for the whole protein pair
	python3 /fsimb/groups/imb-luckgr/imb-luckgr2/projects/AlphaFold/NDD_project_predictions/recycle_precomputed_msas.py \
                                      --prediction_name $prediction_name \
                                      --fasta_path $fasta_path \
                                      --path_to_run $run_path \
                                      --storage_path $store_path \
                                      --precomputed_msas_path $msas_path \
	# document the recycling of the MSAs
	if [ "$(ls -A "${run_dir}${prediction_name}/msas/A/")" ]; then
		echo "${prediction_name}    ${fragment1}" >> recycling_MSA.txt
	else
		if [[ $fragA == "1_0" || $fragA == "1_1" ]]; then
			echo "${prediction_name}    ${fragment1}" >> repeating_MSA.txt
		fi
		if [[ $fragA == "0_0" ]]; then
			echo "${prediction_name}    ${fragment1}" >> first_time_running_MSA.txt
		fi
		if [[ $fragA == "0_1" ]]; then
			echo "${prediction_name}    ${fragment1}" >> malfunctioning_MSA_tracking.txt
		fi
	fi

	if [ "$(ls -A "${run_dir}${prediction_name}/msas/B/")" ]; then
		echo "${prediction_name}    ${fragment2}" >> recycling_MSA.txt
	else
		if [[ $fragB == "1_0" || $fragB == "1_1" ]]; then
			echo "${prediction_name}    ${fragment2}" >> repeating_MSA.txt
		fi
		if [[ $fragB == "0_0" ]]; then
			echo "${prediction_name}    ${fragment2}" >> first_time_running_MSA.txt
		fi
		if [[ $fragB == "0_1" ]]; then
			echo "${prediction_name}    ${fragment2}" >> malfunctioning_MSA_tracking.txt
		fi
	fi
	# AlphaFold Execution
	# However, only run the prediction if it hasn't been already done (check it on the group drive though):
	if [ ! -f "${group_drive_dir}output/runs/run${run_id}/run${prediction_name}/${prediction_name}/ranking_debug.json" ]; then

		# let's change the value in the af_variable.txt file to 1, to indicate that the prediction has started:
		echo "1" > "af_variable_${task_id_here}.txt"

		singularity run --nv --nvccli --writable-tmpfs --bind /home,/fsimb,/media,/mnt,/tmp /mnt/storage/alphafold/v232/alphafold_2.3.2.sif \
                --use_gpu_relax=$GPU_RELAX \
                --model_preset=$MODEL_PRESET \
                --num_multimer_predictions_per_model=1 \
                --max_template_date=2020-05-14 \
                --data_dir=/mnt/storage/alphafold/v232 \
                --bfd_database_path=/mnt/storage/alphafold/v232/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
                --mgnify_database_path=/mnt/storage/alphafold/v232/mgnify/mgy_clusters_2022_05.fa \
                --obsolete_pdbs_path=/mnt/storage/alphafold/v232/pdb_mmcif/obsolete.dat \
                --pdb_seqres_database_path=/mnt/storage/alphafold/v232/pdb_seqres/pdb_seqres.txt \
                --template_mmcif_dir=/mnt/storage/alphafold/v232/pdb_mmcif/mmcif_files \
                --uniprot_database_path=/mnt/storage/alphafold/v232/uniprot/uniprot.fasta \
                --uniref30_database_path=/mnt/storage/alphafold/v232/uniref30/UniRef30_2021_03 \
                --uniref90_database_path=/mnt/storage/alphafold/v232/uniref90/uniref90.fasta \
                --output_dir=$run_dir \
                --fasta_paths=$fasta_file_path \
                --use_precomputed_msas=true \
		> "${OUTPUT_DIR_runs}${output_dir}log_dir/monitoring_dir/std_out_err.log" 2>&1
		# leave use precomputed msas as true just in case

		# change the value in the af_variable.txt file to 0, to indicate that the prediction has finished:
		echo "0" > "af_variable_${task_id_here}.txt"

		# Stage out:
		find "${job_name}/" -maxdepth 6 -mindepth 6 -name "*.pkl" -type  f -exec rm {} \;
		find "${job_name}/" -maxdepth 6 -mindepth 6 -name "unrelaxed_model_*" -type  f -exec rm {} \;
	fi

	find "$run_dir" -maxdepth 4 -mindepth 4 -name "*.sto" -type f -exec rm {} \;
	find "$run_dir" -maxdepth 4 -mindepth 4 -name "*.a3m" -type f -exec rm {} \;

	(
	  flock -x 200  # Exclusive lock on file descriptor 200

	  new_val=1
	  temp_file=$(mktemp) # this creates a unique file
	  awk -v frag1="$fragment1" -v frag2="$fragment2" -v new_val="$new_val" '
	  {
	  if ($1 == frag1 || $1 == frag2) {
	  # Replace the value for either fragment1 or fragment2
	    print $1, new_val
	  } else {
		  # Print unchanged lines
	    print
	  }
	  }' "$finalized_dictionary_file_path" > "$temp_file" && mv "$temp_file" "$finalized_dictionary_file_path"

	) 200>"$finalized_dictionary_file_path.lock"


####################################################################################
else
	if [[ $fragA == "0_1" ]]; then
		echo "${prediction_name}    ${fragment1}" >> malfunctioning_MSA_tracking.txt
	fi
	if [[ $fragB == "0_1" ]]; then
		echo "${prediction_name}    ${fragment2}" >> malfunctioning_MSA_tracking.txt
	fi
	if [[ $fragA == "0_0" ]]; then
		echo "${prediction_name}    ${fragment1}" >> first_time_running_MSA.txt
	fi
	if [[ $fragB == "0_0" ]]; then
		echo "${prediction_name}    ${fragment2}" >> first_time_running_MSA.txt
	fi

	# AlphaFold Execution
	# However, only run the prediction if it hasn't been already done:
	if [ ! -f "${group_drive_dir}output/runs/run${run_id}/run${prediction_name}/${prediction_name}/ranking_debug.json" ]; then

		# let's change the value in the af_variable.txt file to 1, to indicate that the prediction has started:
		echo "1" > "af_variable_${task_id_here}.txt"

		singularity run --nv --nvccli --writable-tmpfs --bind /home,/fsimb,/media,/mnt,/tmp /mnt/storage/alphafold/v232/alphafold_2.3.2.sif \
                --use_gpu_relax=$GPU_RELAX \
                --model_preset=$MODEL_PRESET \
                --num_multimer_predictions_per_model=1 \
                --max_template_date=2020-05-14 \
                --data_dir=/mnt/storage/alphafold/v232 \
                --bfd_database_path=/mnt/storage/alphafold/v232/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
                --mgnify_database_path=/mnt/storage/alphafold/v232/mgnify/mgy_clusters_2022_05.fa \
                --obsolete_pdbs_path=/mnt/storage/alphafold/v232/pdb_mmcif/obsolete.dat \
                --pdb_seqres_database_path=/mnt/storage/alphafold/v232/pdb_seqres/pdb_seqres.txt \
                --template_mmcif_dir=/mnt/storage/alphafold/v232/pdb_mmcif/mmcif_files \
                --uniprot_database_path=/mnt/storage/alphafold/v232/uniprot/uniprot.fasta \
                --uniref30_database_path=/mnt/storage/alphafold/v232/uniref30/UniRef30_2021_03 \
                --uniref90_database_path=/mnt/storage/alphafold/v232/uniref90/uniref90.fasta \
                --output_dir=$run_dir \
                --fasta_paths=$fasta_file_path \
                --use_precomputed_msas=true \
		> "${OUTPUT_DIR_runs}${output_dir}log_dir/monitoring_dir/std_out_err.log" 2>&1
		# leave use precomputed msas as true just in case 

		# change the value in the af_variable.txt file to 0, to indicate that the prediction has finished:
		echo "0" > "af_variable_${task_id_here}.txt"

		# Stage out:
		find "${job_name}/" -maxdepth 6 -mindepth 6 -name "*.pkl" -type  f -exec rm {} \;
		find "${job_name}/" -maxdepth 6 -mindepth 6 -name "unrelaxed_model_*" -type  f -exec rm {} \;
		# this will delete those files in all of the directories though but that's fine

	fi

	find "$run_dir" -maxdepth 4 -mindepth 4 -name "*.sto" -type f -exec rm {} \;
	find "$run_dir" -maxdepth 4 -mindepth 4 -name "*.a3m" -type f -exec rm {} \;
	(
          flock -x 200  # Exclusive lock on file descriptor 200

	  new_val=1
	  temp_file=$(mktemp) # this creates a unique file
	  awk -v frag1="$fragment1" -v frag2="$fragment2" -v new_val="$new_val" '
	  {
	  if ($1 == frag1 || $1 == frag2) {
	  # Replace the value for either fragment1 or fragment2
	    print $1, new_val
	  } else {
		  # Print unchanged lines
	    print
	  }
	  }' "$finalized_dictionary_file_path" > "$temp_file" && mv "$temp_file" "$finalized_dictionary_file_path"
	) 200>"$finalized_dictionary_file_path.lock"

fi

####################################################################################

echo $fragment1 >> all_frag.txt
echo $fragment2 >> all_frag.txt

rm "af_variable_${task_id_here}.txt"


