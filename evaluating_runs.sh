#!/bin/bash

# let's evaluate the first 18 runs that we have for now
python3 model_eval_script.py -run_ids 19,22,78 \
-path_to_run /fsimb/groups/imb-luckgr/imb-luckgr2/projects/AlphaFold/NDD_project_predictions/ndd_project_af_interface_predictions/output/john_predictions/ \
-project_name ndd_project_af_interface_predictions
