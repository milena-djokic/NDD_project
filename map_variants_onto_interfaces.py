# This script maps the 2023 ClinVar, DDD study and Susann's mutation data onto the structures of the interfaces.
# Authors: Milena Djokic, Chop Yan Lee

from pymol import cmd
import os, glob
import pandas as pd
from itertools import cycle
import re

root_path = '/fsimb/groups/imb-luckgr/imb-luckgr2/projects/AlphaFold/NDD_project_predictions/ndd_project_af_interface_predictions/output/john_predictions/'

def parse_prediction_name(prediction_folder):
    """Parse the prediction name of from a prediction folder that has the convention of
    run{run_id}_proteinA_{O/D}_startA_endA.proteinB_{O/D}_startB_endB.

    Args:
        prediction_folder (str): absolute path of the prediction folder

    Returns:
        proteinA (str): uniprot id of protein A
        startA (int): start of fragment A
        endA (int): end of fragment A
        proteinB (str): uniprot id of protein B
        startB (int): start of fragment B
        endB (int): end of fragment B
    """
    prediction_name = os.path.split(prediction_folder)[-1]
    prediction_name = '_'.join(prediction_name.split('_')[1:])
    regionA, regionB = prediction_name.split('.')
    proteinA = regionA.split('_')[0]
    startA, endA = regionA.split('_')[2:]
    proteinB = regionB.split('_')[0]
    startB, endB = regionB.split('_')[2:]
    return proteinA, int(startA), int(endA), proteinB, int(startB), int(endB)


def load_fragment_fasta_sequence(file_path, uniprot_id, fragment_start, fragment_end):
    sequence = None
    load_next = False
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):  # Header line
                if uniprot_id in line and str(fragment_start) + '_' + str(fragment_end) in line:
                    load_next = True
                else:
                    load_next = False
            elif load_next:
                sequence = line
                break
    return sequence


def retrieve_variant_info(uniprot_id,fragment_start,fragment_end):
    """Retrieve protein-level variants from 2023 ClinVar, DDD study and Susann Schweiger mutation data.
       Group it based on NDD association and pathogenicity.
    Args:
        uniprot_id (str): Uniprot ID of a protein
        fragment_start (int): Start of the predicted fragment
        fragment_end (int): End of the predicted fragment

    Return:
        variant_dict (dict): Dictio
    """
    mut = pd.read_csv('mut.tsv', header=0, sep='\t')
    # we will turn conflicting mutations to uncertain, for simplicity reasons:
    mut.replace({'conflicting':'uncertain'}, inplace=True)
    mut['mut_id'] = [a+str(b)+c for a,b,c in zip(mut['ref'],mut['pos'],mut['alt'])]
    protein_mut = mut[mut['uniprot_id'] == uniprot_id] # in case we ever need this for anything
    fragment_mut = protein_mut[(protein_mut['pos'] >= fragment_start) & (protein_mut['pos'] <= fragment_end)]

    # now we will split them up into 6 dictionaries:
    pathogenic_ndd = fragment_mut[(fragment_mut['pathogenicity'] == 'pathogenic') & (fragment_mut['ndd_phe'] == 1)]
    pathogenic_non_ndd = fragment_mut[(fragment_mut['pathogenicity'] == 'pathogenic') & (fragment_mut['ndd_phe'] == 0)]
    uncertain_ndd = fragment_mut[(fragment_mut['pathogenicity'] == 'uncertain') & (fragment_mut['ndd_phe'] == 1)]
    uncertain_non_ndd = fragment_mut[(fragment_mut['pathogenicity'] == 'uncertain') & (fragment_mut['ndd_phe'] == 0)]
    benign_ndd = fragment_mut[(fragment_mut['pathogenicity'] == 'benign') & (fragment_mut['ndd_phe'] == 1)]
    benign_non_ndd = fragment_mut[(fragment_mut['pathogenicity'] == 'benign') & (fragment_mut['ndd_phe'] == 0)]
    # I have to think about the situations when one or all of them are empty and how to deal with this
    all_mut_dict = {mutation:[fragment_start, fragment_end, protein_pos, protein_pos-fragment_start] for mutation, protein_pos in zip(fragment_mut['mut_id'],fragment_mut['pos'])}
    pathogenic_ndd_dict = {mutation:[fragment_start, fragment_end, protein_pos, protein_pos-fragment_start] for mutation, protein_pos in zip(pathogenic_ndd['mut_id'],pathogenic_ndd['pos'])}
    pathogenic_non_ndd_dict = {mutation:[fragment_start, fragment_end, protein_pos, protein_pos-fragment_start] for mutation, protein_pos in zip(pathogenic_non_ndd['mut_id'],pathogenic_non_ndd['pos'])}
    uncertain_ndd_dict = {mutation:[fragment_start, fragment_end, protein_pos, protein_pos-fragment_start] for mutation, protein_pos in zip(uncertain_ndd['mut_id'],uncertain_ndd['pos'])}
    uncertain_non_ndd_dict = {mutation:[fragment_start, fragment_end, protein_pos, protein_pos-fragment_start] for mutation, protein_pos in zip(uncertain_non_ndd['mut_id'],uncertain_non_ndd['pos'])}
    benign_ndd_dict = {mutation:[fragment_start, fragment_end, protein_pos, protein_pos-fragment_start] for mutation, protein_pos in zip(benign_ndd['mut_id'],benign_ndd['pos'])}
    benign_non_ndd_dict = {mutation:[fragment_start, fragment_end, protein_pos, protein_pos-fragment_start] for mutation, protein_pos in zip(benign_non_ndd['mut_id'],benign_non_ndd['pos'])}

    variant_dict = {'all_variants':all_mut_dict,'pathogenic_ndd':pathogenic_ndd_dict, 'pathogenic_non_ndd':pathogenic_non_ndd_dict, 'uncertain_ndd':uncertain_ndd_dict, 'uncertain_non_ndd':uncertain_non_ndd_dict, 'benign_ndd':benign_ndd_dict, 'benign_non_ndd':benign_non_ndd_dict}

    return variant_dict

def annotate_variant_info(prediction_folder):
    """Use retrieve_variant_info to retrieve all protein-level variants from all three mutation datasets for both of the fragments.
       Launch already existing .pse files (generated by the script model_eval_script.py), to map the mutation info on the interfaces.
       Mutations are saved as possibly 6 objects (colored based on the category) - pathogenic_ndd, pathogenic_non_ndd, uncertain_ndd, uncertain_non_ndd, benign_ndd and benign_non_ndd.
       The file is saved as _variant_annotated for viewing.

       Args:
           prediction_folder (str): absolute path of the prediction folder
    """

    # load the _contacts.pse model for variant annotation
    processed_pse_files = glob.glob(f'{prediction_folder}/model_eval/ranked_0_*_variant_annotated.pse')

    for processed_file in processed_pse_files:
        if os.path.exists(processed_file):
            return

    for file in glob.glob(f'{prediction_folder}/model_eval/ranked_0_*.pse'):

        cmd.load(file)
        wt_object = f'{prediction_folder}'.split('/')[-1]
        fragA_start = wt_object.split('.')[0].split('_')[-2] # we will need this for later
        fragB_start = wt_object.split('.')[1].split('_')[-2]

        cmd.set_color('darkred', (98,12,12))
        cmd.set_color('lightred', (255,92,92))
        cmd.set_color('darkyellow', (219,161,0))
        cmd.set_color('lightyellow', (255,223,132))
        cmd.set_color('darkblue', (18,13,156))
        cmd.set_color('lightblue', (186,183,255))
        color_dictionary = {'pathogenic_ndd':'darkred','pathogenic_non_ndd':'lightred','uncertain_ndd':'darkyellow','uncertain_non_ndd':'lightyellow','benign_ndd':'darkblue','benign_non_ndd':'lightblue'}

        no_mutations = 0
        for chain_id, variant_dict in zip(['A','B'],[variant_dictA,variant_dictB]):

            # if there are no mutations in one of the chains, just continue to the next chain
            if not variant_dict.get('all_variants'):
               no_mutations += 1
               continue

            if chain_id == 'A':
                frag_start = fragA_start
                color_backbone = 'green'
            if chain_id == 'B':
                frag_start = fragB_start
                color_backbone = 'magenta'

            for variant_type in list(variant_dict.keys())[1:]:

                # if we don't have any variants of this particular type, skip it
                if len(variant_dict[variant_type]) == 0:
                    continue

                # if we have more than one variant of the same type on a particular residue, label them as such:
                position_dictionary = {}
                for position in sorted(set([value[-1] for value in variant_dict[variant_type].values()])):
                    position_dictionary[position] = [key for key, value in variant_dict[variant_type].items() if value[-1] == position]

                stick_color = color_dictionary.get(variant_type)
                for variant in variant_dict[variant_type]:
                    variant_pos = variant[1:-1]
                    variant_pos_ini = int(variant_pos) - int(frag_start)
                    variant_pos = variant_pos_ini + 1

                    wt_aa = variant[0]
                    mut_aa = variant[-1]
                    AA_dict = {'Phe':'F','Leu':'L','Ser':'S','Tyr':'Y','Cys':'C','Trp':'W','Pro':'P','His':'H','Gln':'Q','Arg':'R','Ile':'I','Met':'M','Thr':'T','Asn':'N','Lys':'K','Val':'V','Ala':'A','Asp':'D','Glu':'E','Gly':'G','Ter':'*'}
                    AA_dict = {v: k for k, v in AA_dict.items()}
                    wt_aa_long = AA_dict.get(wt_aa)
                    mut_aa_long = AA_dict.get(mut_aa)
                    aa_label = f'{wt_aa_long} -> {mut_aa_long}'

                    cmd.select(name=f'{variant_type}_var', selection=f'chain {chain_id} and resi {variant_pos}')

                    # create an object that will contain all the amino acids of this type
                    if not cmd.get_names("objects").count(f"{variant_type}"):
                        # If the object doesn't exist, create it
                        cmd.create(f"{variant_type}", f"{variant_type}_var", zoom=0)

                        # color the sticks for this residue
                        cmd.color(stick_color, f'{variant_type}') # first we basically color everything this color, this is because otherwise the C atoms in the sticks remain green
                        cmd.set('cartoon_color', color_backbone, f'chain {chain_id}') # then we leave the backbone green or magenta!
                        cmd.color('atomic', 'not elem C') # then we color atoms as they should be
                        cmd.set('stick_color', stick_color, f'{variant_type} and not elem H and not elem N and not elem O and not elem S') # then we color the stick, except the elements
                        cmd.label(f'({variant_type} and chain {chain_id} and resi {variant_pos} and name CA)', f'"{aa_label}"')

                        #cmd.label(f'chain {chain_id} and {variant_type} and resi {variant_pos} and name CA', f'"{aa_label}"')

                    else:
                        # If it exists, merge the selection
                        cmd.create(f"{variant_type}", f"{variant_type} or {variant_type}_var", zoom=0)

                        # color the sticks for this residue
                        cmd.color(stick_color, f'{variant_type}') # first we basically color everything this color, this is because otherwise the C atoms in the sticks remain green
                        cmd.set('cartoon_color', color_backbone, f'chain {chain_id}') # then we leave the backbone green or magenta!
                        cmd.color('atomic', 'not elem C') # then we color atoms as they should be
                        cmd.set('stick_color', stick_color, f'{variant_type} and not elem H and not elem N and not elem O and not elem S') # then we color the stick, except the elements

                        # we only have to do this here:
                        if len(position_dictionary.get(variant_pos_ini)) > 1:
                            wt_aa_long = AA_dict.get(position_dictionary.get(variant_pos_ini)[0][0])
                            mut_1_aa_long = AA_dict.get(position_dictionary.get(variant_pos_ini)[0][-1])
                            remaining_mut_list = []
                            for aa in position_dictionary.get(variant_pos_ini)[1:]:
                                remaining_mut_list.append(AA_dict.get(aa[-1]))
                            remaining_mut = ' | '.join(remaining_mut_list)
                            aa_label = f'{wt_aa_long} -> {mut_1_aa_long} | {remaining_mut}'
                            cmd.label(f'({variant_type} and chain {chain_id} and resi {variant_pos} and name CA)', f'"{aa_label}"')

                        else:
                            cmd.label(f'({variant_type} and chain {chain_id} and resi {variant_pos} and name CA)', f'"{aa_label}"')

                    cmd.show_as("sticks", f"{variant_type}")
                    for selection in cmd.get_names('selections'):
                        if selection != 'residues_less_5A':
                            cmd.delete(selection)

        # we want the sequence indices to reflect the fragment start and ends
        cmd.alter("chain A", f"resi=str(int(resi) + {fragA_start} - 1)")
        cmd.alter("chain B", f"resi=str(int(resi) + {fragB_start} - 1)")
        cmd.rebuild()

        # reorder the objects:
        obj_order = ['pathogenic_ndd','pathogenic_non_ndd','uncertain_ndd','uncertain_non_ndd','benign_ndd','benign_non_ndd']
        obj_order = [i for i in obj_order if i in cmd.get_names('objects')]
        cmd.order(' '.join(obj_order))

        cmd.show("labels", "all")

        if no_mutations == 2:
            file_name = f'{file[:-4]}_no_variant_annotated.pse'
        else:
            file_name = f'{file[:-4]}_variant_annotated.pse'
        cmd.save(file_name)


if __name__ == '__main__':

    #runs = [i+1 for i in range(33)]
    runs = [19, 22, 78]
    for run in runs:
        run_path = os.path.join(root_path, f'run{run}')
        dir_names = [name for name in os.listdir(run_path) if os.path.isdir(os.path.join(run_path, name))]

        smiley_faces = ['(｡♥‿♥｡)','\(^ _ ^)/','(∩˃o˂∩)♡','♥‿♥','(　◠ ◡ ◠　)','♪(๑ᴖ◡ᴖ๑)♪','(๑^ں^๑)','(*＾v＾*)','☝( ◠‿◠ )☝','(︶ω︶)','（＾＿－）']
        for smiley, prediction in zip(cycle(smiley_faces), dir_names):
            print(f'{smiley} Yaaaay we got another one down! {prediction}')

            proteinA, startA, endA, proteinB, startB, endB = parse_prediction_name(prediction)

            # fetch the variant data for both fragments
            variant_dictA = retrieve_variant_info(proteinA,startA,endA)
            variant_dictB = retrieve_variant_info(proteinB,startB,endB)

            pymol_path = os.path.join(run_path, prediction)
            annotate_variant_info(pymol_path)

            # find the fasta file of the prediction so that we can also compare the sequence with the mutation and if it's not the same, raise a flag:
            fasta_path = '/fsimb/groups/imb-luckgr/imb-luckgr2/projects/AlphaFold/NDD_project_predictions/ndd_project_af_interface_predictions/input/john_predictions/'
            fasta_path = os.path.join(fasta_path, re.findall(r'\d+', prediction.split('_')[0])[0], prediction[3:]+'.fasta')
            fragA_seq = load_fragment_fasta_sequence(fasta_path, proteinA, startA, endA)
            fragB_seq = load_fragment_fasta_sequence(fasta_path, proteinB, startB, endB)

            # compare to the fasta sequence and raise a flag if there is a mismatch!
            #for key, value in zip(variant_dictA['all_variants'].keys(), variant_dictA['all_variants'].values()):
            #    ref_aa = key[0]
            #    pos_aa = value[-1]
            #    print(ref_aa, pos_aa)
            #    print(fragA_seq, fragA_seq[int(pos_aa)])
            #    if fragA_seq[int(pos_aa)] != ref_aa:
            #        print(f'Mutation does not match the sequence! Prediction: {prediction}, fragment: {proteinA}_{startA}_{endA}, mutation: {key},{value}')

            #for key, value in zip(variant_dictB['all_variants'].keys(), variant_dictB['all_variants'].values()):
            #    ref_aa = key[0]
            #    pos_aa = value[-1]
            #    print(ref_aa, pos_aa)
            #    print(fragB_seq, fragB_seq[int(pos_aa)])
            #    if fragB_seq[int(pos_aa)] != ref_aa:
            #        print(f'Mutation does not match the sequence! Prediction: {prediction}, fragment: {proteinB}_{startB}_{endB}, mutation: {key},{value}')
















