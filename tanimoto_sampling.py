# Copyright (c) 2022-2026 Samuel Plamondon

# Performance monitoring:

import time
start = time.time()

# Imports:

import csv
import pandas as pd
from random import shuffle
from rdkit import Chem
from rdkit import DataStructs

# Functions:

def import_rct_list(in_file):
    
    """
    in_file: filename of reactant list (incl. extension)
    returns: tuple containing list of SMILES strings (smiles_list) and list of
        compound ID strings (id_list)
    """
    
    df = pd.read_csv(in_file, sep = ',', dtype = str)
    smiles_list = df[smiles_head]
    id_list = df[id_head]
    
    return smiles_list, id_list

def randomize_list_order(smiles_list, id_list):
    
    """
    smiles_list: list of SMILES strings
    id_list: list of compound IDs (strings)
    returns: a tuple containing shuffled SMILES and ID lists (shuff_smiles,
        shuff_id)
    """
    
    shuff_smiles = []
    shuff_id = []
    indices = list(range(len(smiles_list)))
    shuffle(indices)
    
    for i in indices:
        
        shuff_smiles.append(smiles_list[i])
        shuff_id.append(id_list[i])
    
    return shuff_smiles, shuff_id

def tanimoto_trim(shuff_smiles, shuff_id, tani_threshold):
    
    """
    shuff_smiles: list of SMILES strings
    shuff_id: list of compound IDs (strings)
    tani_threshold: float, upper limit for allowed Tanimoto similarity
    returns: a tuple containing trimmed SMILES and ID lists (kept_smiles,
        kept_id)
    """
    
    kept_smiles = [shuff_smiles[0]]
    kept_id = [shuff_id[0]]
    shuff_mol = [Chem.MolFromSmiles(smiles) for smiles in shuff_smiles]
    fingers = [Chem.RDKFingerprint(shuff_mol[0])]
    
    for i in range(1, len(shuff_mol)):
        
        finger = Chem.RDKFingerprint(shuff_mol[i])
        too_similar = False
        
        for f in fingers:
            
            if DataStructs.FingerprintSimilarity(finger, f) > tani_threshold:
                
                too_similar = True
                break
        
        if not too_similar:
            
            kept_smiles.append(shuff_smiles[i])
            kept_id.append(shuff_id[i])
            fingers.append(finger)
    
    return kept_smiles, kept_id

def tani_threshold_opt(shuff_smiles, shuff_id, target_size):
    
    """
    shuff_smiles: list of SMILES strings
    shuff_id: list of compound IDs (strings)
    target_size: int, desired number of molecules in trimmed list
    returns: a tuple containing (a) trimmed SMILES and ID lists of strings;
        and (b) a float corresponding to the final percent error of the list
        length vs. target length (kept_smiles, kept_id, percent_error)
    """
    
    target_proportion = target_size/len(shuff_smiles)
    min_threshold = 0.0
    max_threshold = 1.0
    err_threshold = 10
    percent_error = err_threshold + 1 # Arbitrary value to allow the following
                                      # "while" loop to initiate.
    
    while (percent_error > err_threshold and
           (max_threshold - min_threshold) > 0.001):
        
        tani_threshold = (max_threshold + min_threshold)/2
        
        print("Trying Tanimoto threshold of " + str(tani_threshold) + "... ",
              end = "", flush = True)
        
        kept_smiles, kept_id = tanimoto_trim(
            shuff_smiles, shuff_id, tani_threshold)
        test_proportion = len(kept_smiles)/len(shuff_smiles)
        percent_error = 100*abs(test_proportion/target_proportion - 1)
        
        print(str(round(percent_error, 1)) + "%" + " off target.")
        
        if test_proportion > target_proportion:
            max_threshold = tani_threshold
        else:
            min_threshold = tani_threshold
        
    print("")
    
    return kept_smiles, kept_id, percent_error

def write_rct_file(kept_smiles, kept_id, out_file):
    
    """
    kept_smiles: list of SMILES strings
    kept_id: list of compound IDs (strings)
    out_file: output file name
    Writes the SMILES and ID strings to a CSV file, returns None.
    """
    
    output_smiles = [smiles_head]
    output_id = [id_head]
    
    output_smiles += kept_smiles
    output_id += kept_id
    
    with open(out_file, "a", newline = '') as output_file:
        
        writer = csv.writer(output_file)
        writer.writerows(zip(output_smiles, output_id))
    
    output_file.close

# Inputs:

# Reads an input file, "input_tanimoto_sampling.txt", with the following lines:

# Inputs
# Header of SMILES column
# Header of ID column
# Number of lists (n)
# Target list size
# Input filename (with extension) 1
# Input filename (with extension) 2
# ...
# Input filename (with extension) n
# Output filename (with extension) 1
# Output filename (with extension) 2
# ...
# Output filename (with extension) n

inputs = pd.read_csv('input_tanimoto_sampling.txt', sep = ',', dtype = str)
inputs = inputs['Inputs']

smiles_head = inputs[0]
id_head = inputs[1]
num_lists = int(inputs[2])
target_size = int(inputs[3])

in_file = []
out_file = []

for x in range(num_lists):
    
    in_file.append(inputs[4 + x])
    out_file.append(inputs[4 + x + num_lists])

# Main program:

for i in range(num_lists):
    
    smiles_list, id_list = import_rct_list(in_file[i])
    
    start_num = len(smiles_list)
    shuff_smiles, shuff_id = randomize_list_order(smiles_list, id_list)
    
    del(smiles_list)
    del(id_list)
    
    kept_smiles, kept_id, percent_error = tani_threshold_opt(
        shuff_smiles, shuff_id, target_size)
    
    end_num = len(kept_smiles)
    
    write_rct_file(kept_smiles, kept_id, out_file[i])
    
    print("List #" + str(i + 1) + ": Out of " + str(start_num) +
          " molecules, " + str(end_num) + " were kept. The target was " + 
          str(target_size) + ".")
    print("")

# Completion of performance monitoring:

end = time.time()
tottime = (end - start)/60
print(f"Runtime (min): {tottime}")
