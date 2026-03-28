# This script automates the organization of MSAs that are created by AlphaFold so as to enable AlphaFold to reuse the MSAs that it has generated for subsequent predictions that involve the same fragments. Important: This script only works if the prediction file name is standardized, i.e. run118_O75381_O_105_224.P50542_O_1_322.
# Author: Milena Djokic and Chop Yan Lee
import os, json, argparse, sys, glob, shutil, string

def parse_args(args=sys.argv[1:]):
    """Parse arguments"""
    parser = argparse.ArgumentParser()
    # the current prediction for which we are running this script:
    parser.add_argument(
        "--prediction_name", type=str, help="Prediction name", dest="curr_pred_name"
    )
    # path to where we keep the fasta files
    parser.add_argument(
        "--fasta_path", type=str, help="Fasta path", dest="fasta_path"
    )
    # path to where we run our predictions - in my case it's gonna be the scratch space
    parser.add_argument(
        "-run_path", "--path_to_run", type=str, help="Run path", dest="run_path"
    )
    # path to where we store the data resulting from the predictions
    parser.add_argument(
        "-store_path", "--storage_path", type=str, help="Storage path", dest="store_path"
    )
    # path to where we store the computed MSAs - central directory
    parser.add_argument(
        "-computed_msas_path", "--precomputed_msas_path", type=str, help="Msas path", dest="computed_msas_path"
    )
    return parser.parse_args(args)

def make_msas_dir_chain_id_map(run_path,fasta_path,computed_msas_path): # so here the prediction path is actually the run path!
    """Make a prediction directory of yet to be predicted fasta file. In this prediction directory, make also an msas directory with chain_id_map.json in the msas directory.
    
    Args:
        run_path (str): the absolute path of the prediction folder
        fasta_path (str): the absolute path of the folder where we store the fasta file
        computed_msas_path (str): the absolute path where the computed MSAs are stored (the central MSAs directory)
    """
    msas_path = f'{run_path}msas/'
    path_to_prediction, prediction_name = os.path.split(run_path)

    # make the msas folder in the prediction folder if it has not been created from the previous organizing of msas folders yet
    if not os.path.exists(msas_path):
        os.makedirs(msas_path)

    # parse prediction fasta file
    descriptions = []
    sequences = []
    index = -1
    with open(fasta_path, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
    for line in lines:
        if line.startswith('>'):
            index += 1
            descriptions.append(line[1:])
            sequences.append('')
            continue
        sequences[index] += line

    chain_id_map = {}
    for chain_id, description, sequence in zip(string.ascii_uppercase,descriptions,sequences):
        chain_id_map[chain_id] = {"description":description,"sequence":sequence}

    # make chain_id_map.json file if not exist
    if not os.path.exists(f'{msas_path}chain_id_map.json'):
        with open(f'{msas_path}chain_id_map.json','w') as f:
            json.dump(chain_id_map,f,indent=4)

    # with the chain_id_map, check if msas of prediction exist in central directory
    for chain_id in chain_id_map:
        chain_name = chain_id_map.get(chain_id).get('description')
        uniprot_id = chain_name.split('_')[0]
        central_dir = f'{computed_msas_path}{uniprot_id}/'
        # if msas already copied before, skip this chain_id
        if os.path.exists(f'{msas_path}{chain_id}/mgnify_hits.sto'):
            continue
        # if msas files already exist in central dir, move it to the current prediction folder
        if os.path.exists(f'{central_dir}{chain_name}/mgnify_hits.sto'):
            if not os.path.exists(f'{msas_path}{chain_id}'):
                os.makedirs(f'{msas_path}{chain_id}',exist_ok=True)
            for msas_file in glob.glob(f'{central_dir}{chain_name}/*'):
                shutil.copy(msas_file,f'{msas_path}{chain_id}/')

def move_computed_msas_to_central(store_path,computed_msas_path):
    """Move the msas computed in the previous AlphaFold prediction to a central directory where future predictions involving the same fragment can be looked up. Check if the prediction folder has been predicted, if so check if central directory has all the msa files, if all the msas already exist in the central directory, remove the msas from the prediction folder to not duplicate msas.
    
    Args:
        store_path (str): the absolute path to where we store the prediction folder
        computed_msas_path (str): the absolute path where the computed MSAs are stored (the central MSAs directory)
    """
    msas_path = f'{store_path}/msas/'
    path_to_prediction, prediction_name = os.path.split(store_path)
    # print(f'Organizing msas of {prediction_name} in central directory...')

    with open(f'{msas_path}chain_id_map.json','r') as f:
        chain_id_map = json.load(f)
    for chain_id in chain_id_map:
        chain_name = chain_id_map.get(chain_id).get('description')
        uniprot_id = chain_name.split('_')[0]
        central_dir = f'{computed_msas_path}{uniprot_id}/'
        # check if msas do not exist in central directory, move there
        if not os.path.exists(f'{central_dir}{chain_name}/'):
        	os.makedirs(f'{central_dir}{chain_name}',exist_ok=True)
        # check every msa file if they exist in the central directory, if it doesn't exist, move there, else, remove to save storage
        # there might be cases where the MSA files have already been moved to another prediction folder, and then this will just skip
        for msas_file in glob.glob(f'{msas_path}{chain_id}/*'):
            _, msas = os.path.split(msas_file)
            if not os.path.exists(f'{central_dir}{chain_name}/{msas}'):
            	shutil.copy(msas_file,f'{central_dir}{chain_name}/')
            else:
                pass

if __name__ == "__main__":
    
    args = parse_args()
    curr_pred_name = vars(args)["curr_pred_name"]
    fasta_path = vars(args)["fasta_path"]
    run_path = vars(args)["run_path"]
    store_path = vars(args)["store_path"]
    computed_msas_path = vars(args)["computed_msas_path"]
    fasta_files = [file for file in glob.glob(fasta_path + '/*.fasta')]
    

    #print(f"Child process PID: {os.getpid()}")
    predicted = []
    not_predicted = [curr_pred_name]
    for file in fasta_files:
    	prediction_name = file.split("/")[-1]
    	prediction_name = '.'.join(prediction_name.split('.')[:2])
    	run_subdirectory = f'run{prediction_name}'
    	# we should not be checking if the prediction is done but if the MSA files are there or not!
    	cond1=os.path.exists(f'{store_path}{run_subdirectory}/{prediction_name}/msas/A/mgnify_hits.sto')
    	cond2=os.path.exists(f'{store_path}{run_subdirectory}/{prediction_name}/msas/B/mgnify_hits.sto')
    	if cond1 or cond2: #predicted or not, MSA files are there!
        	predicted.append(prediction_name)
    	else:
        	pass
    # for those already predicted, move their msas to the central directory and delete them from the prediction folder
    for pred in predicted:
    	# now we have to make the store path and computed_msas_path prediction specific:
    	store_path_pred = f'{store_path}run{pred}/{pred}/'
    	move_computed_msas_to_central(store_path_pred,computed_msas_path)

    # for those not yet predicted, check if their msas exist in central directory and if so move to the prediction folder
    if curr_pred_name != "None":
    	for pred in not_predicted:
    		run_path_pred = f'{run_path}run{pred}/{pred}/'
    		fasta_path_pred = f'{fasta_path}{pred}.fasta'
    		make_msas_dir_chain_id_map(run_path_pred,fasta_path_pred,computed_msas_path)
    else:
    	pass

#fasta_path="${INPUT_DIR_runs}${run_id}/" # input dir from the scratch directory
#run_path="${OUTPUT_DIR_runs}run${run_id}/" # output dir in the scratch directory 
#store_path="${group_drive_dir}output/runs/run${run_id}/" # stored predictions
#msas_path="${group_drive_dir}output/computed_msas/" # central directory for recycling msas
## and these dirs are not run/prediction specific but rather for the whole protein pair
#python3 /fsimb/groups/imb-luckgr/imb-luckgr2/projects/AlphaFold/NDD_project_predictions/recycle_precomputed_msas.py \
#                                      --prediction_name $prediction_name \
#                                      --fasta_path $fasta_path \
#                                      --path_to_run $run_path \
#                                      --storage_path $store_path \
#                                      --precomputed_msas_path $msas_path
