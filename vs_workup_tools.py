##Searches a directory tree for xml files and strips docking scores. Aggregates them in a single file and ranks them according to score in ascending order

import xml.etree.ElementTree as ET
import os
import glob
import pandas as pd
import numpy as np
import time
import seaborn as sns
import matplotlib.pyplot as plt
import csv
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
import pathlib
import os
import shutil
import glob
import time
import json
from functools import partial
import subprocess
import csv
from ast import literal_eval
import shlex
import re
import collections
from rdkit import Chem
import datetime;
import math

def compound_indexer(target_dir, subjob_input_db_ext, run_dir):
    csv_columns = ['index_id', 'compound_name']
    with open(f'{run_dir}/compound_index_key.csv', 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
        writer.writeheader()
        file_ext = subjob_input_db_ext
        # print(target_dir)
        for file in os.listdir(target_dir):
            if file.endswith(file_ext):
                file_path = (os.path.join(target_dir, file))
                # print(file_path)
                suppl = Chem.SDMolSupplier(file_path)
                for m in suppl:
                    compound_name = m.GetProp('_Name')
                    index_id = os.path.basename(file_path)
                    # print(index_id)
                    # print(compound_name)
                    d = dict();
                    d['index_id'] = index_id
                    d['compound_name'] = compound_name
                    # print(d)
                    writer.writerow(d)
                    # time.sleep(1)
def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def compound_indexer_v2(temp_dir, stage_dir, subjob_input_db_ext, num_batches):
    sdf_path = []
    for file in os.listdir(temp_dir):
        if file.endswith(subjob_input_db_ext):
            file_path = (os.path.join(temp_dir, file))
            suppl = Chem.SDMolSupplier(file_path)
            for m in suppl:
                canonical_id = m.GetProp('_Name')
                index_id = os.path.basename(file_path)
                sdf_path.append(tuple([file_path, canonical_id]))
    batches_list = list(enumerate(sdf_path))
    # print(sdf_path)
    # len(batches_list)
    batch_dir_num = int(len(batches_list)/num_batches)
    # print(batch_dir_num)
    chunked_lis = chunks(sdf_path, (batch_dir_num+1))
    chunked_batch = dict(enumerate(chunked_lis))
    # len(chunked_batch)
    canonical_id = []
    sn_id = []
    input_file_path = []
    output_batch_path = []
    input_file_ext = []

    for i in range(len(chunked_batch)):
        for targets in chunked_batch[i]:
            compound_id = (pathlib.Path(targets[0])).stem
            input_path = os.path.join(stage_dir, r'input/batch{:04d}/'.format(i+1), compound_id + subjob_input_db_ext)
    #         print(input_path)
            output_batch_path2 = os.path.join(stage_dir, r'output/batch{:04d}/'.format(i+1))
    #         print(targets[0], destination, targets[1])
            canonical_id.append(targets[1])
            sn_id.append(compound_id)
            input_file_path.append(input_path)
            output_batch_path.append(output_batch_path2)
            input_file_ext.append(subjob_input_db_ext)
    result5 = zip(sn_id, output_batch_path, input_file_path, input_file_ext)
    canonical_params = list(result5)
    result6 = dict(zip(canonical_id, canonical_params))

    # print(result6)
    jsonfilepath = os.path.join(stage_dir + 'dictionary.json')
    with open(jsonfilepath, 'w', encoding='utf-8') as jsonf:
        jsonf.write(json.dumps(result6, indent=4))
    return jsonfilepath

def dump_autodockgpu_scores_tocsv(run_dir):
    rootdir = run_dir
    names = [os.path.basename("*.csv")]
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            #print os.path.join(subdir, file)
            filepath = subdir + os.sep + file
            if filepath.endswith(".xml"):
                # print(filepath)
                tree = ET.parse(filepath)
                root = tree.getroot()
                with open(filepath+'dockingscore.csv','w') as outfile:
                    for thing in root.findall(".//*[@cluster_rank='1']"):
                        things_baby = thing.get('mean_binding_energy')
                        if things_baby is not None:
                                outfile.write("{}\n".format(things_baby))

def import_compound_index_to_df(run_dir):
    index_key = os.path.join(run_dir, 'compound_index_key.csv')
    # print(f'the location of the index key is {index_key}')    
    file_df = pd.read_csv(index_key)
    file_df['index_id'] = file_df['index_id'].map(lambda x: x.rstrip('.sdf'))
    file_df = file_df.sort_values('index_id',ascending=True)
    file_df.rename(columns={'index_id':'names'},inplace=True)
    # print(file_df)
    return(file_df)


def compile_autodock_gpu_scores_in_df(run_dir):
    rootdir = run_dir
    df = pd.DataFrame()
    #for file_ in all_files:
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
        #print os.path.join(subdir, file)
            filepath = subdir + os.sep + file
            if filepath.endswith("xmldockingscore.csv"):
                file2_df = pd.read_csv(filepath,sep=';', infer_datetime_format=True,header=None )
                file2_df['names'] = file
                df = df.append(file2_df)
        

    df.rename(columns={0:'Docking_Score'},inplace=True)
    df['Rank'] = df['Docking_Score'].rank(method='first')
    df = df.set_index('Rank')
    df = df.sort_values('names',ascending=True)
    df['names'] = df['names'].map(lambda x: x.rstrip('.dlg.xmldockingscore.csv'))
    compiled_raw_scores = os.path.join(run_dir, 'raw_docking_score_ranks.csv')
    df.to_csv(compiled_raw_scores) 
    # print(df)
    return df

def merge_compound_id_name_and_plot(run_dir, compound_index_df, autodockscores_df):

    plt.style.use('ggplot')
    df_merged = pd.merge(autodockscores_df, compound_index_df, on='names', how='outer')
    # print(df_merged)
    df_rank_sorted = df_merged.sort_values('Docking_Score',ascending=True)
    # print(df_rank_sorted)

    docking_output = os.path.join(run_dir, 'docking_score_dist.pdf')
    sns.histplot(data=df_rank_sorted, x='Docking_Score')
    plt.savefig('docking_score_dist.pdf', dpi=300)

    scores_output = os.path.join(run_dir, 'docking_score_ranks.csv')
    df_rank_sorted.to_csv(scores_output)





def extract_top_poses_from_stage(target_dir):
    rootdir = target_dir
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            #print os.path.join(subdir, file)
            filepath = Path(subdir + os.sep + file)
#             print(filepath)
            if filepath.suffix == '.xml':
                # print(filepath)
                target_run = find_topscoring_run(filepath)
                pose_file = filepath.with_suffix('.dlg')
                # print(pose_file)
                extract_top_autodock_pose(pose_file, target_run)
                

def find_topscoring_run(path):
    tree = ET.parse(path)
    root = tree.getroot()
    for data in root.findall(".//*[@cluster_rank='1']"):
        target_run = data.get('run')
        if target_run is not None:
            return target_run

def extract_top_autodock_pose(file_name, target_run):
    """
    Gets the first pdbqt structure of .dlg to create a new pdbqt
    file with the same file name.
    """

    file = open(file_name, "r")
    lines = file.readlines()
    target_dir = str(file_name.parent)
    
    
    starting_line = 0
    ending_line = 0
    # print(rf'Run:   {target_run} / 10')
    for line in lines:
        if line.startswith(rf'Run:   {target_run} / 10'):
            starting_line = lines.index(line)
            break
            
    if target_run != 10:
        # print('the target was not 10')
        for line in lines[starting_line:]:
            if line.startswith(rf'Run:   {(int(target_run) + 1)} / 10'):
                ending_line = lines.index(line)
                break
        pdbqt_name_prefix = file_name.stem + ".pdbqt"
        pdbqt_name = os.path.join(target_dir, pdbqt_name_prefix)
        pdbqt_content = lines[(starting_line + 4):(ending_line - 9)]
        stripped_pdbqt_content = [line.rstrip('\n') for line in pdbqt_content]

        
    if target_run == 10:
        # print('the target was 10')
        for line in lines[starting_line:]:
            if line.startswith('    CLUSTERING HISTOGRAM'):
                ending_line = lines.index(line)
                break
        pdbqt_name_prefix = file_name.stem + ".pdbqt"
        pdbqt_name = os.path.join(target_dir, pdbqt_name_prefix)
        pdbqt_content = lines[(starting_line + 4):(ending_line - 5)]
        stripped_pdbqt_content = [line.rstrip('\n') for line in pdbqt_content]
        
    clean_pdbqt_content = []
    
    for line in stripped_pdbqt_content:
        cleaned_line = line[max(line.find('D'), 8):]
        if not cleaned_line.startswith('USER'):
            clean_pdbqt_content.append(cleaned_line)
        
    with open(pdbqt_name, 'w') as f:
        for line_item in clean_pdbqt_content:
            f.write("%s\n" % line_item)


    # df = pd.DataFrame()
    # #for file_ in all_files:
    # for subdir, dirs, files in os.walk(rootdir):
    #     for file in files:
    #         #print os.path.join(subdir, file)
    #         filepath = subdir + os.sep + file
    #         if filepath.endswith("_docking_scores.csv"):
    #             file_df = pd.read_csv(filepath,sep=';', infer_datetime_format=True,header=None )
    #             file_df['names'] = file
    #             df = df.append(file_df)
        
    # df.rename(columns={0:'Docking_Score'},inplace=True)
    # df['Rank'] = df['Docking_Score'].rank(method='first')
    # df = df.set_index('Rank')
    # df = df.sort_values('Docking_Score',ascending=True)
    # df['names'] = df['names'].map(lambda x: x.rstrip('.xmldockingscore.csv'))
    # print(df)

    # df.to_csv(f'{run_dir}/{run_name}_docking_scores.csv')


# #### may need to download numpy in order to get this part to work.
# import pandas as pd
# import numpy as np

# df = pd.read_csv('/Volumes/home/sweetnothings/various_tests/testscores.csv', sep = ',')
# def func(a):
#     if "actives" in a.lower():
#         return "yes"
#     elif "decoys" in a.lower():
#         return "no"
#     else:
#         return "Other"

# df["column2_type"] = df.names.apply(lambda x: func(x))

# for column2_type, data in df.groupby('column2_type'):
#     data.to_csv("{}.csv".format(column2_type))
# print (df)
# df.to_csv('.csv')
