import os
import time
import shutil
# import click
import math
import typing
import glob 
from typing import Set
from cleanup_utils import cleanup_scripts_dir, cleanup_dir, cleanup_input_dir
from rich.console import Console
import csv
from pydantic import BaseModel
import json
import config
import molchunk_config
from files_folders_utils2 import resubmit_failed_jobs, check_failed_job_status, input_vs_output_count, trackJobs, check_job_status, make_batch_submit_script, make_submit_script, master_chunk_depositor, setup_run_dir, next_path, makefolders, list_files, make_stage_dir, make_temp_dir, submit_sge_scripts
from molecule_functions import split_molecule_database, MolCount
from vs_workup_tools import extract_top_poses_from_stage, compound_indexer, dump_autodockgpu_scores_tocsv, import_compound_index_to_df, compile_autodock_gpu_scores_in_df, merge_compound_id_name_and_plot
import shutil
import pathlib
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
import datetime;
import math
from datetime import timedelta
from timeit import time
import re
import xml.etree.ElementTree as ET
from ast import literal_eval
import shlex
import collections
from shutil import copytree,copy2
from openeye import oechem
from tqdm.contrib import tenumerate
from tqdm import tqdm
import pyarrow.dataset as ds
import pyarrow.parquet as pq

console = Console()

class Stage(BaseModel):
    name: str
    desc: str
    node_type: str
    num_cores: int
    setup_fn: typing.Any = None
    move_files_fn: typing.Any = None
    cleanup_fn: typing.Any = None
    exit_fn: typing.Any = None
    exe_cmd: str
    input_db_ext: str
    output_db_ext: str
    dependency_header: str


class Pipeline(BaseModel):
    name: str
    working_dir: str
    scratch: str
    file_prefix: str
    stages: typing.List[Stage]
    molecule_db: str
    molecule_db_ext: str
    output_db_ext: str
    


class Job(BaseModel):
    pass

def stopwatch(method):
    def timed(*args, **kw):
        ts = time.perf_counter()
        result = method(*args, **kw)
        te = time.perf_counter()
        duration = timedelta(seconds=te - ts)
        print(f"{method.__name__}: {duration}")
        return result
    return timed

@stopwatch
def prepare_db_fn(input_file, num_parts, stage_dir, scratch, ext, SDF_name):
    temp_prefix_path = os.path.join(scratch, SDF_name)
    print(f'Writing SDF chunks to {scratch}')
    outbase = temp_prefix_path

    moldb = oechem.OEMolDatabase(input_file)
    molcount = moldb.NumMols()

    chunksize, lft = divmod(molcount, num_parts)
    if lft != 0:
        chunksize += 1
    chunk, count = 1, 0

    ofs = create_output_stream(outbase, ext, chunk)
    for idx in range(moldb.GetMaxMolIdx()):
        count += 1
        if count > chunksize:
            if chunk == lft:
                chunksize -= 1

            ofs.close()
            chunk, count = chunk + 1, 1
            ofs = create_output_stream(outbase, ext, chunk)

        moldb.WriteMolecule(ofs, idx)


def create_output_stream(outbase, ext, chunk):
    newname = outbase + ('_%05d' % chunk) + ext
    ofs = oechem.oemolostream()
    if not ofs.open(newname):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % newname)
    return ofs

@stopwatch
def move_files_between_stages(stage_dir, stage_count, run_dir, temp_dir, input_db_ext):
    dict_path = os.path.join(run_dir, f'{pathlib.Path(stage_dir).parent.stem}_dictionary.json') 

    with open(dict_path) as json_file: 
        data = json.load(json_file)

    previous_stage_dir = os.path.join(run_dir, f'stage_{stage_count}')

    console.print(f'The previous stage directory was: {previous_stage_dir}', style="bold underline")

    for l in data['canonical_ids']:      
        try:
            assigned_batch_key = (pathlib.Path(previous_stage_dir).stem+'_batch')
            previous_output_file_path = os.path.join(previous_stage_dir, 'output', data['canonical_ids'][l][assigned_batch_key], data['canonical_ids'][l]['job_assigned_id']+input_db_ext)
            temp_file_location = os.path.join(temp_dir, data['canonical_ids'][l]['job_assigned_id']+input_db_ext)
            shutil.move(previous_output_file_path, temp_file_location)
        except FileNotFoundError:
            print("Wrong file or file path")

    return previous_stage_dir

    # previous_stage_dir = os.path.join(run_dir, f'stage_{stage_count}')
    # target_dir = current_stage_temp_dir
    # source_dir = os.path.join(previous_stage_dir, 'output')
    # l_of_directories = glob.glob(os.path.join(source_dir, "*", ""))
    # # print(f'The list of directories in previous output are: {l_of_directories}')
    # for d in range(len(l_of_directories)):
    #     batch_dir = l_of_directories[d]
    #     # print(f'The batch_dir is: {batch_dir}')
    #     file_names = glob.glob(os.path.join(batch_dir, "*" + subjob_input_db_ext))
    #     # print(f'The file names are: {file_names}')
    #     for file_name in file_names:
    #         shutil.copy(os.path.join(file_name), target_dir)
    # return source_dir
@stopwatch
def setup_stage(run_dir, stage_name, batch_size):
    print(f'The run directory path is: {run_dir}')
    stage_dir = make_stage_dir(run_dir, stage_name)
    print(f'The stage directory is: {stage_dir}')
    time.sleep(2)
    
    input_batch_directories = [r'input/batch{:04d}'.format(i) for i in range(1, int(batch_size) + 1)]
    makefolders(stage_dir, input_batch_directories)

    scripts_batch_directories = [r'scripts/batch{:04d}/'.format(i) for i in range(1, int(batch_size) + 1)]
    makefolders(stage_dir, scripts_batch_directories)
    
    output_batch_directories = [r'output/batch{:04d}'.format(i) for i in range(1, int(batch_size) + 1)]
    makefolders(stage_dir, output_batch_directories)
    
    temp_dir = make_temp_dir(stage_dir)
    
    stage_dirs = {
        'input_batch_dir': os.path.join(stage_dir + '/input/'),
        'scripts_batch_dirs': scripts_batch_directories,
        'temp_dir': temp_dir,
        'stage_dir': stage_dir,
        'output_batch_dir': output_batch_directories,
        'input_batch_dirs': input_batch_directories
    }
    return stage_dirs
@stopwatch
def copy2_verbose(src, dst):
    print('Copying {0}'.format(src))
    copy2(src,dst)

@stopwatch
def cleanup_stage_dir(run_dir, stage_dir, temp_dir, scripts_dir, input_dir, batch_size):
    cleanup_scripts_dir(run_dir, stage_dir, scripts_dir, batch_size)
    cleanup_input_dir(run_dir, stage_dir, input_dir, batch_size)
    cleanup_dir(temp_dir)

@stopwatch
def prep_docking_results(run_dir, stage_dir, output_db_ext):
    dict_path = os.path.join(run_dir, f'{pathlib.Path(stage_dir).parent.stem}_dictionary.json') 
    with open(dict_path) as json_file: 
        data = json.load(json_file)

    for l in data['canonical_ids']:
        assigned_batch_key = (pathlib.Path(stage_dir).stem+'_batch')
        job_assigned_id = data['canonical_ids'][l]['job_assigned_id']
        output_dir = os.path.join(stage_dir, 'output')
        xml_path = os.path.join(output_dir, data['canonical_ids'][l][assigned_batch_key], f'{job_assigned_id}{output_db_ext}.xml')
        dlg_path = os.path.join(output_dir, data['canonical_ids'][l][assigned_batch_key], f'{job_assigned_id}{output_db_ext}.dlg')
        # print(dlg_path)
        try:
            tree = ET.parse(xml_path)
            root = tree.getroot()
            docking_score = [x.get('mean_binding_energy') for x in root.findall(".//*[@cluster_rank='1']")]
            top_docking_run = [y.get('run') for y in root.findall(".//*[@cluster_rank='1']")]
            # print(docking_score)
            # print(top_docking_run)
            data['canonical_ids'][l].update({'mean_binding_energy':docking_score, 'top_docking_run':top_docking_run, 'combined_docking_poses':dlg_path})
        except FileNotFoundError:
            print(f'There was an error with prep_docking_results')
            print(f'The target xml_path was: {xml_path}')
            print("Wrong file or file path")

    with open(dict_path, 'w', encoding='utf-8') as jsonf:
        jsonf.write(json.dumps(data, indent=4))


@stopwatch
def extract_top_autodock_pose(run_dir, stage_dir):     
    
    dict_path = os.path.join(run_dir, f'{pathlib.Path(stage_dir).parent.stem}_dictionary.json') 
    with open(dict_path) as json_file: 
        data = json.load(json_file)

    for l in data['canonical_ids']:
        try:
            file_name = pathlib.Path(data['canonical_ids'][l]['combined_docking_poses'])
            target_run_list = data['canonical_ids'][l]['top_docking_run']
            target_run = target_run_list[0]
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
                pre_pdbqt_name_prefix = os.path.splitext(file_name.stem)
                pdbqt_name_prefix = pre_pdbqt_name_prefix[0] + ".pdbqt"
                pdbqt_name = os.path.join(target_dir, pdbqt_name_prefix)
                pdbqt_content = lines[(starting_line + 4):(ending_line - 9)]
                stripped_pdbqt_content = [line.rstrip('\n') for line in pdbqt_content]

                
            if target_run == 10:
                # print('the target was 10')
                for line in lines[starting_line:]:
                    if line.startswith('    CLUSTERING HISTOGRAM'):
                        ending_line = lines.index(line)
                        break
                pre_pdbqt_name_prefix = os.path.splitext(file_name.stem)
                pdbqt_name_prefix = pre_pdbqt_name_prefix[0] + ".pdbqt"
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
        except KeyError:
            print(f'There was a key error for extract_top_autodock_pose')
            print(f'The target dlg path was : {file_name}')
            pass

@stopwatch
def workup_docking_results(run_dir, stage_dir, output_db_ext):
    prep_docking_results(run_dir, stage_dir, output_db_ext)
    extract_top_autodock_pose(run_dir, stage_dir)

@stopwatch
def index_3d_pose(run_dir, stage_dir, output_db_ext):
    dict_path = os.path.join(run_dir, f'{pathlib.Path(stage_dir).parent.stem}_dictionary.json') 
    with open(dict_path) as json_file: 
        data = json.load(json_file)

        for l in data['canonical_ids']:
            try:
                output_dir = os.path.join(stage_dir, 'output')
                assigned_batch_key = (pathlib.Path(stage_dir).stem+'_batch')
                job_assigned_id = data['canonical_ids'][l]['job_assigned_id']
                pdb_path = os.path.join(output_dir, data['canonical_ids'][l][assigned_batch_key], f'{job_assigned_id}{output_db_ext}')
                data['canonical_ids'][l].update({'final_docking_pose':pdb_path})
            except KeyError:
                pass

    with open(dict_path, 'w', encoding='utf-8') as jsonf:
        jsonf.write(json.dumps(data, indent=4))

@stopwatch
def indexer(run_dir, stage_dir, scratch, input_db_ext, batch_size, num_parts):
    console.print(f'The assigned batch size for this stage is {batch_size}', style="bold underline")
    run_id = pathlib.Path(stage_dir).parent.stem
    sdf_path = [os.path.join(scratch, run_id+'_{:05d}'.format(i)+input_db_ext) for i in range(1,num_parts+1)]
    # print(sdf_path)
    # sdf_path = glob.glob(os.path.join(temp_dir, f'*{input_db_ext}'))
    sdf_per_batch = math.ceil(len(sdf_path)/batch_size)
    total_canon_id = len(sdf_path)
    canon_id_per_batch = math.ceil(total_canon_id/batch_size)


    Dict = {}
    Dict['canonical_ids'] = {}
    for count, sdf in enumerate(sdf_path):
        with open(sdf, 'rb') as reader:

            suppl = Chem.ForwardSDMolSupplier(reader)
            canonical_id = ' '.join(map(str, [m.GetProp('_Name') for m in suppl]))

            index_id = pathlib.Path(sdf).stem
                
            #     batch_num = math.ceil(count/canon_id_per_batch)
            assigned_batch = 'batch{:04d}'.format(math.ceil(count/canon_id_per_batch))
            assigned_batch_key = (pathlib.Path(stage_dir).stem+'_batch')
            # print(assigned_batch)
            # print(assigned_batch_key)
            if math.ceil(count/canon_id_per_batch) == 0:
                assigned_batch = 'batch{:04d}'.format(1)
            else:
                assigned_batch = 'batch{:04d}'.format(math.ceil(count/canon_id_per_batch))
            Dict['canonical_ids'][canonical_id] = {'job_assigned_id':index_id, assigned_batch_key:assigned_batch}



    jsonfilepath = os.path.join(run_dir, f'{pathlib.Path(stage_dir).parent.stem}_dictionary.json')

    with open(jsonfilepath, 'w', encoding='utf-8') as jsonf:
        jsonf.write(json.dumps(Dict, indent=4))
    return jsonfilepath

@stopwatch
def dataset_indexer(run_dir, stage_dir, input_db_ext, batch_size, num_parts, fragments, program_exe):
    console.print(f'The assigned batch size for this stage is {batch_size}', style="bold underline")
    run_id = pathlib.Path(stage_dir).parent.stem
    # sdf_path = [os.path.join(scratch, run_id+'_{:05d}'.format(i)+input_db_ext) for i in range(1,num_parts+1)]
    # print(sdf_path)
    # sdf_path = glob.glob(os.path.join(temp_dir, f'*{input_db_ext}'))
    # sdf_per_batch = math.ceil(len(sdf_path)/batch_size)
    # total_canon_id = len(sdf_path)
    canon_id_per_batch = math.ceil(num_parts/batch_size)


    Dict = {}
    Dict['canonical_ids'] = {}
    for count, df_frag in tenumerate(fragments, start=0, total=num_parts):
        # with open(sdf, 'rb') as reader:

        # suppl = Chem.ForwardSDMolSupplier(reader)
        # canonical_id = ' '.join(map(str, [m.GetProp('_Name') for m in suppl]))
        canonical_id = run_id + '_' + str(count)
        index_id = canonical_id
            
        #     batch_num = math.ceil(count/canon_id_per_batch)
        assigned_batch = 'batch{:04d}'.format(math.ceil(count/canon_id_per_batch))
        assigned_batch_key = (pathlib.Path(stage_dir).stem+'_batch')
        # print(assigned_batch)
        # print(assigned_batch_key)
        if math.ceil(count/canon_id_per_batch) == 0:
            assigned_batch = 'batch{:04d}'.format(1)
        else:
            assigned_batch = 'batch{:04d}'.format(math.ceil(count/canon_id_per_batch))
        

        input_file_path = os.path.join(stage_dir, 'input', assigned_batch, index_id+input_db_ext)

        # table = df_frag.to_table()

        # pq.write_table(table, input_file_path)

        output_file_path = os.path.join(stage_dir, 'output', assigned_batch_key, index_id+input_db_ext)

        exe_cmd = program_exe.format(input_file_path=pathlib.Path(input_file_path).parent, output_file_path=pathlib.Path(output_file_path).parent, prefix=run_id)

        exe_cmd_key = pathlib.Path(stage_dir).stem

        Dict['canonical_ids'][canonical_id] = {'job_assigned_id':index_id, assigned_batch_key:assigned_batch, exe_cmd_key:exe_cmd}

        # data['canonical_ids'][l].update({exe_cmd_key:exe_cmd})

    jsonfilepath = os.path.join(run_dir, f'{pathlib.Path(stage_dir).parent.stem}_dictionary.json')

    with open(jsonfilepath, 'w', encoding='utf-8') as jsonf:
        jsonf.write(json.dumps(Dict, indent=4))
    return jsonfilepath


@stopwatch
def batch_indexer(dict_path, stage_dir, temp_dir, input_db_ext, batch_size, num_parts):
    console.print(f'The assigned batch size for this stage is {batch_size}', style="bold underline")
    run_id = pathlib.Path(stage_dir).parent.stem
    sdf_path = [os.path.join(temp_dir, run_id+'_{:05d}'.format(i)+input_db_ext) for i in range(1,num_parts+1)]
    # sdf_path = glob.glob(os.path.join(temp_dir, f'*{input_db_ext}'))
    sdf_per_batch = math.ceil(len(sdf_path)/batch_size)
    total_canon_id = len(sdf_path)
    canon_id_per_batch = math.ceil(total_canon_id/batch_size)

    with open(dict_path) as json_file: 
        data = json.load(json_file)

    for count, l in enumerate(data['canonical_ids']):
            #     batch_num = math.ceil(count/canon_id_per_batch)
            assigned_batch = 'batch{:04d}'.format(math.ceil(count/canon_id_per_batch))
            assigned_batch_key = (pathlib.Path(stage_dir).stem+'_batch')
            # print(assigned_batch)
            # print(assigned_batch_key)
            if math.ceil(count/canon_id_per_batch) == 0:
                assigned_batch = 'batch{:04d}'.format(1)
            else:
                assigned_batch = 'batch{:04d}'.format(math.ceil(count/canon_id_per_batch))
            data['canonical_ids'][l].update({assigned_batch_key:assigned_batch})



    # jsonfilepath = os.path.join(pathlib.Path(stage_dir).parent.stem, f'{pathlib.Path(stage_dir).parent.stem}_dictionary.json')

    with open(dict_path, 'w', encoding='utf-8') as jsonf:
        jsonf.write(json.dumps(data, indent=4))
    return dict_path

@stopwatch
def shuttler(dict_path, temp_dir, stage_dir, program_exe, input_db_ext, output_db_ext):

    with open(dict_path) as json_file: 
        data = json.load(json_file)

    for l in data['canonical_ids']:
        try:
            assigned_batch_key = (pathlib.Path(stage_dir).stem+'_batch')
            temp_file_location = os.path.join(temp_dir, data['canonical_ids'][l]['job_assigned_id']+input_db_ext)
            input_file_path = os.path.join(stage_dir, 'input', data['canonical_ids'][l][assigned_batch_key], data['canonical_ids'][l]['job_assigned_id']+input_db_ext)
            shutil.move(temp_file_location, input_file_path)
            output_file_path = os.path.join(stage_dir, 'output', data['canonical_ids'][l][assigned_batch_key], data['canonical_ids'][l]['job_assigned_id']+output_db_ext)
            exe_cmd = program_exe.format(input_file_path=input_file_path, output_file_path=output_file_path)
            exe_cmd_key = pathlib.Path(stage_dir).stem
            data['canonical_ids'][l].update({exe_cmd_key:exe_cmd})
        except FileNotFoundError:
            console.print("Warning: Shuttler unabale to move files, Wrong file or file path due to indexing issue.", )
            console.print(f'Shuttler tried to move: {temp_file_location} to this location: {input_file_path}')
        


    with open(dict_path, 'w', encoding='utf-8') as jsonf:
        jsonf.write(json.dumps(data, indent=4))
    return dict_path

@stopwatch
def scratch_shuttler(dict_path, scratch, stage_dir, program_exe, input_db_ext, output_db_ext):

    with open(dict_path) as json_file: 
        data = json.load(json_file)

    for l in data['canonical_ids']:
        try:
            assigned_batch_key = (pathlib.Path(stage_dir).stem+'_batch')
            temp_file_location = os.path.join(scratch, data['canonical_ids'][l]['job_assigned_id']+input_db_ext)
            input_file_path = os.path.join(stage_dir, 'input', data['canonical_ids'][l][assigned_batch_key], data['canonical_ids'][l]['job_assigned_id']+input_db_ext)
            shutil.move(temp_file_location, input_file_path)
            output_file_path = os.path.join(stage_dir, 'output', data['canonical_ids'][l][assigned_batch_key], data['canonical_ids'][l]['job_assigned_id']+output_db_ext)
            exe_cmd = program_exe.format(input_file_path=input_file_path, output_file_path=output_file_path)
            exe_cmd_key = pathlib.Path(stage_dir).stem
            data['canonical_ids'][l].update({exe_cmd_key:exe_cmd})
        except FileNotFoundError:
            console.print("Warning: Shuttler unabale to move files, Wrong file or file path due to indexing issue.", )
            console.print(f'Shuttler tried to move: {temp_file_location} to this location: {input_file_path}')
        


    with open(dict_path, 'w', encoding='utf-8') as jsonf:
        jsonf.write(json.dumps(data, indent=4))
    return dict_path

# @stopwatch
# def script_maker(dict_path, stage_dir):
#     with open(dict_path) as json_file: 
#         data = json.load(json_file)

#     scripts_dir = os.path.join(stage_dir, 'scripts') 
#     for l in data['canonical_ids']:
#         try:
#             assigned_batch_key = (pathlib.Path(stage_dir).stem+'_batch')
#             input_prefix = data['canonical_ids'][l]['job_assigned_id']
#             job_filename = f'run_{input_prefix}.sh'
#             script_file_path = os.path.join(scripts_dir, data['canonical_ids'][l][assigned_batch_key], job_filename)
#             exe_cmd_key = pathlib.Path(stage_dir).stem
#             exe_cmd = data['canonical_ids'][l][exe_cmd_key]
#             with open(script_file_path, 'w') as job_f:
#                 shell_script = '''\
#                 {exe_cmd}
#                 '''.format(exe_cmd=exe_cmd)
#                 job_f.write(shell_script)
#         except KeyError:
#             print('keyerror')
#             pass

@stopwatch
def next_path(path_pattern):
    """
    Finds the next free path in an sequentially named list of files
    e.g. path_pattern = 'file-%s.txt':
    file-1.txt
    file-2.txt
    file-3.txt
    Runs in log(n) time where n is the number of existing files in sequence
    """
    i = 1

    # First do an exponential search
    while os.path.exists(path_pattern % i):
        i = i * 2

    # Result lies somewhere in the interval (i/2..i]
    # We call this interval (a..b] and narrow it down until a + 1 = b
    a, b = (i // 2, i)
    while a + 1 < b:
        c = (a + b) // 2 # interval midpoint
        a, b = (c, b) if os.path.exists(path_pattern % c) else (a, c)

    return path_pattern % b



def make_dataset_submit_script(stage_dir, scripts_batch_directories, cpu_or_gpu, job_dependency, input_db_ext, run_dir, program_exe, input_dataset):
    csv_columns = ['run', 'stage', 'batch', 'inputbatchpath', 'numberoffilesin', 'sgecommand', 'sgejobnumber']
    # # csv_file = print(f'{print(os.path.basename(run_dir))}.csv')
    csv_file = os.path.join(run_dir, os.path.basename(stage_dir) + '_submission_metadata.csv')
    with open(csv_file, 'w') as csvfile:
        
        writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
        writer.writeheader()
        for batch_directories in scripts_batch_directories:
            batch_directory_path = os.path.join(stage_dir, batch_directories)
            # read_files = glob.glob(os.path.join(batch_directory_path, "*.sh"))
            file_path = os.path.join(batch_directory_path, 'submit_jobs.sh')
            run_id = pathlib.Path(stage_dir).parent.stem
            batch_number = os.path.basename(os.path.normpath(batch_directory_path))
            input_batch_path = os.path.join(input_dataset, batch_number)
            output_batch_path = os.path.join(stage_dir, 'output', batch_number)

            exe_cmd = program_exe.format(input_file_path=input_batch_path, output_file_path=output_batch_path, prefix=run_id)
            #print(batch_directory_path)
            with open(file_path, "w") as outfile:
                outfile.write(job_dependency)
                shell_script = '''\
                {exe_cmd}
                '''.format(exe_cmd=exe_cmd)
                outfile.write(shell_script)

            #print(file_path)
            


            job_submit_info = submit_sge_scripts(cpu_or_gpu, input_batch_path, file_path, input_db_ext, stage_dir, run_dir)

            # ee = add_job_parameters(file_path, batch_directory_path, cpu_or_gpu, stage_dir, job_dependency, input_db_ext, run_dir)
        
            writer.writerow(job_submit_info)





@stopwatch
def run():
    STAGES = [
        Stage(
            name='stage_1',
            desc='Enumerate SMILES into stereoisomers and generate fingerprints',
            node_type='cpu',
            num_cores=150,
            setup_fn=None,
            move_files_fn=None,
            cleanup_fn=cleanup_stage_dir,
            exit_fn=None,
            exe_cmd='python /cbica/home/grahamth/molbeam/molbeam/app2.py {input_file_path} {output_file_path}',
            input_db_ext='.parquet',
            output_db_ext='.parquet',
            dependency_header='''\

            /bin/echo "Loading miniconda 3.8 and starting ray_test release enviroment..."
            source /cbica/software/external/python/miniconda/3.8/etc/profile.d/conda.sh
            module load python/miniconda/3.8
            conda activate ray_test
            export OE_LICENSE=/cbica/home/grahamth/oe_license.txt
            export PATH=/cbica/home/grahamth/openeyeapp_2020_2/openeye/bin:$PATH

            /bin/echo "Checking OE_License location..."
            echo $OE_LICENSE

            /bin/echo "Checking conda enviroment..."
            which conda

            /bin/echo "Unsetting PYTHONPATH..."
            unset PYTHONPATH

            /bin/echo "Checking python..."
            which python

            /bin/echo "PYTHON env details"
            env | grep PYTHON
            ''',
        ),
        # Stage(
        #     name='stage_2',
        #     desc='Conversion of PDB files to PDBQT files.',
        #     node_type='cpu',
        #     num_cores=300,
        #     setup_fn=None,
        #     move_files_fn=move_files_between_stages,
        #     cleanup_fn=cleanup_stage_dir,
        #     exit_fn=None, 
        #     exe_cmd='/cbica/home/grahamth/autodock/ADFRsuite_x86_64Linux_1.0/buildtest/bin/prepare_ligand -v -l {input_file_path} -o {output_file_path}',
        #     input_db_ext='.pdb',
        #     output_db_ext='.pdbqt',
        #     dependency_header='''\
        #     pwd
        #     ''',
        # ),
        # Stage(
        #     name='stage_3',
        #     desc='Docking of PDBQT files to receptor grid.',
        #     node_type='gpu',
        #     num_cores=200,
        #     setup_fn=None,
        #     move_files_fn=move_files_between_stages,
        #     cleanup_fn=cleanup_stage_dir,
        #     exit_fn=workup_docking_results,
        #     exe_cmd='''\
        #     /cbica/home/grahamth/autodockgpucuda10_2/AutoDock-GPU/bin/autodock_gpu_128wi \
        #             -ffile /cbica/home/grahamth/autodock/dopamine/dude_d3/pocket2_fixer_moreatoms/rigidReceptor.maps.fld \
        #             -lfile {input_file_path} \
        #             -resnam {output_file_path} \
        #             -nrun 10
        #     ''',
        #     input_db_ext='.pdbqt',
        #     output_db_ext='.dlg',
        #     dependency_header='''\
        #     module unload cuda
        #     module load cuda/10.2
        #     CUDA_VISIBLE_DEVICES=$(get_CUDA_VISIBLE_DEVICES) || exit
        #     export CUDA_VISIBLE_DEVICES
        #     /bin/echo "CUDA_VISIBLE_DEVICES"
        #     echo $CUDA_VISIBLE_DEVICES
        #     /bin/echo "nvcc --version"
        #     nvcc --version
        #     /bin/echo "nvidia-smi"
        #     nvidia-smi
        #     ''',
        # ),
        # Stage(
        #     name='stage_4',
        #     desc='Conversion of top docking pose PDBQT to PDB files',
        #     node_type='cpu',
        #     num_cores=300,
        #     setup_fn=None,
        #     move_files_fn=move_files_between_stages,
        #     cleanup_fn=cleanup_stage_dir,
        #     exit_fn=index_3d_pose,
        #     exe_cmd='obabel {input_file_path} -O {output_file_path}',
        #     input_db_ext='.pdbqt',
        #     output_db_ext='.pdb',
        #     dependency_header='''\
        #     export PATH=$PATH:/cbica/home/grahamth/openbabel3/bin
        #     export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cbica/home/grahamth/openbabel3/lib
        #     ''',
        # ),
    ]
    PIPELINE = Pipeline(
        name='er_real_0_to_10_fp_app2_take2_tau/',
        working_dir='/cbica/home/grahamth/sn_production/fp_search',
        scratch='/scratch/grahamth',
        file_prefix='er_real_0_to_10_app2_prod_test/',
        stages=STAGES,
        molecule_db="/cbica/home/grahamth/sweetnothings/various_tests/er_real_0_to_10/stage_1/output",
        molecule_db_ext='parquet',
        output_db_ext='.parquet',           

    )

    print('                                                                                     ')
    print('                                                                                     ')
    console.print(f'Initiating {len(STAGES)}-stage workflow.', style="bold underline")
    print('                                                                                     ')
    print('                                                                                     ')
    # time.sleep(2)
    print('=====================================================================================')
    print('                                                                                     ')
    print(f'Stage 1 = {STAGES[0].desc}                                                          ')
    print('                                                                                     ')
    # # time.sleep(2)
    # print(f'Stage 2 = {STAGES[1].desc}                                                          ')
    # # time.sleep(2)
    # print('                                                                                     ')
    # print(f'Stage 3 = {STAGES[2].desc}                                                          ')
    # # time.sleep(2)
    # print('                                                                                     ')
    # print(f'Stage 4 = {STAGES[3].desc}                                                          ')
    # # time.sleep(2)
    print(f'Current working directory = {PIPELINE.working_dir}                                  ')
    print(f'Parent multicompound SDF = {PIPELINE.molecule_db}                                   ')
    print('                                                                                     ')
    print('                                                                                     ')
    print('                                                                                     ')
    print('=====================================================================================')
    print('=====================================================================================')
    # time.sleep(2)

    print('                                                                                     ')
    print('                                                                                     ')
    # num_parts = MolCount(PIPELINE.molecule_db)
    # print(f'Parent multicompound SDF contains = {num_parts} molecules.                          ')
    
    # dataset = ds.dataset(PIPELINE.molecule_db, format=PIPELINE.molecule_db_ext)
    # fragments = [file for file in dataset.get_fragments()]
    input_dataset = PIPELINE.molecule_db
    print(input_dataset)
    print(glob.glob(f'{input_dataset}/batch*'))
    num_batches = len(glob.glob(f'{input_dataset}/batch*'))

    run_dir = setup_run_dir(PIPELINE.working_dir, PIPELINE.name)
    for count, stage in enumerate(PIPELINE.stages):
        print('                                                                                     ')
        print('                                                                                     ')      
        print(f'Starting {stage.name}.                                                              ')
        print(f'{stage.desc}.                                                                       ')
        print(f'Creating directories, moving files and generating submit scripts.                   ')
        dirs = setup_stage(run_dir, stage.name, num_batches)
        if stage.setup_fn != None:
            # setup_fn = stage.dataset_indexer
            print('                                                                                     ')
            print('                                                                                     ')
            print('=====================================================================================')
            print('=====================================================================================')
            print('                                                                                     ')
            print('                                                                                     ')
            # scratch_temp = os.path.join(PIPELINE.scratch, PIPELINE.name+'_tmp_%s')
            # unique_scratch_temp = next_path(scratch_temp)
            # os.mkdir(unique_scratch_temp)
            # print(f'Setting up the temp directory at: {unique_scratch_temp}.')
            print(f'Setting up {stage.name}. Importing parquet database and splitting into batches')
            # setup_fn(PIPELINE.molecule_db, 
            #     num_parts, 
            #     dirs['stage_dir'], 
            #     unique_scratch_temp, 
            #     PIPELINE.molecule_db_ext, 
            #     PIPELINE.name)
            batch_size = stage.num_cores
            print(f'Indexing target structures.')
            dict_path = dataset_indexer(run_dir, dirs['stage_dir'], stage.input_db_ext, batch_size, num_parts, fragments, stage.exe_cmd)
            # print(f'Moving indexed structures from the scratch to stage_1 in the local working directory: {PIPELINE.working_dir}')
            # scratch_shuttler(dict_path, unique_scratch_temp, dirs['stage_dir'], stage.exe_cmd, stage.input_db_ext, stage.output_db_ext)
            # shutil.rmtree(unique_scratch_temp)
            # print('The scratch temp directory was deleted.')
        print('                                                                                     ')
        print('                                                                                     ')
        # print(f'The previous stage was : Stage_{count}')
        
        # temp_dir_path = os.path.join(dirs['temp_dir'] + '/')

        print('                                                                                     ')
        print('                                                                                     ')
        print('=====================================================================================')
        print('=====================================================================================')
        print('                                                                                     ')
        print('                                                                                     ')

        # if stage.move_files_fn != None:
        #     move_files_between_stages_fn = stage.move_files_fn
        #     previous_stage_output_dir = move_files_between_stages_fn(dirs['stage_dir'], count, 
        #         run_dir, 
        #         temp_dir_path, 
        #         stage.input_db_ext)
        #     stage_dir = dirs['stage_dir']
        #     dict_path = os.path.join(run_dir, f'{pathlib.Path(stage_dir).parent.stem}_dictionary.json')
        #     batch_size = stage.num_cores
        #     batch_indexer(dict_path, dirs['stage_dir'], dirs['temp_dir'], stage.input_db_ext, batch_size, num_parts)
        #     shuttler(dict_path, dirs['temp_dir'], dirs['stage_dir'], stage.exe_cmd, stage.input_db_ext, stage.output_db_ext)
        #     cleanup_dir(previous_stage_output_dir)
        #     print('The previous stage output dir was cleaned up.')
        # batch_size = stage.num_cores
        # molecules_per_batch = math.ceil(num_parts / batch_size)
        
        

        # script_maker(dict_path, dirs['stage_dir'])


        if stage.node_type == 'cpu':
            print('CPU resources have selected for this stage.')

            make_dataset_submit_script(
                dirs['stage_dir'], 
                dirs['scripts_batch_dirs'], 
                molchunk_config.SGE_CPU_PARAM, 
                stage.dependency_header, 
                stage.input_db_ext, 
                run_dir,
                stage.exe_cmd,
                input_dataset
            )
        # elif(stage.node_type == 'gpu'):
        #     print('GPU resources have been slected for this stage.')
        #     make_batch_submit_script(
        #         dirs['stage_dir'], 
        #         dirs['scripts_batch_dirs'], 
        #         config.SGE_GPU_PARAM, 
        #         stage.dependency_header, 
        #         stage.input_db_ext, 
        #         run_dir
        #     )
        else:
            print('Invalid resource request: Please choose cpu or gpu.')

        print('=====================================================================================')
        print('=====================================================================================')
        print('                                                                                     ')
        print('Jobs successfully submitted to SGE. Tracking job status now before starting next stage.')
        print('                                                                                     ')
        submittied_job_list = check_job_status(run_dir, dirs['stage_dir'])
        print(f'The list of SGE assigned jobids are: {submittied_job_list}                          ')
        
        # print(dirs['stage_dir'])
        # print(dirs['input_batch_dirs'])
        # print(dirs['output_batch_dir'])
        # trackJobs(submittied_job_list, 
        #     dirs['stage_dir'], 
        #     num_batches, 
        #     dirs['input_batch_dirs'], 
        #     dirs['output_batch_dir'], 
        #     stage.input_db_ext, 
        #     stage.output_db_ext, 
        #     waittime=15)

        print('                                                                                     ')
        print('                                                                                     ')
        print(f'Checking for job failures and attempting to resubmit.                               ')
        # resubmit_failed_jobs(dirs['output_batch_dir'], 
        #     dirs['stage_dir'], 
        #     run_dir, 
        #     stage.input_db_ext, 
        #     stage.output_db_ext)
        # failed_job_status = check_failed_job_status(run_dir, dirs['stage_dir'])
        # print(failed_job_status)
        # if failed_job_status == 0:
        #     print('There are no failed jobs.')

        # else:
        #     resubmitted_job_list = check_failed_job_status(run_dir, dirs['stage_dir'])
        #     trackJobs(resubmitted_job_list, 
        #         dirs['stage_dir'], 
        #         num_batches, 
        #         dirs['input_batch_dirs'], 
        #         dirs['output_batch_dir'], 
        #         stage.input_db_ext, 
        #         stage.output_db_ext, 
        #         waittime=15)
        # if stage.exit_fn != None:
        #     exit_stage_fn = stage.exit_fn
        #     exit_cleanup_fn = stage.cleanup_fn
        #     print('                                                                                     ')
        #     print('                                                                                     ')
        #     print('=====================================================================================')
        #     print('=====================================================================================')
        #     print('                                                                                     ')
        #     print('                                                                                     ')
        #     print(f'Carrying out exit activities for {stage.name}.                                      ')
        #     exit_stage_fn(run_dir, dirs['stage_dir'], stage.output_db_ext)
        #     input_dir = os.path.join(dirs['stage_dir'], 'input')
        #     scripts_dir = os.path.join(dirs['stage_dir'], 'scripts')
        #     exit_cleanup_fn(run_dir, dirs['stage_dir'], dirs['temp_dir'], scripts_dir, input_dir, batch_size)
        #     print('                                                                                     ')
        #     print('Exit functions complete.                                                             ')
        #     print('=====================================================================================')
        #     print('=====================================================================================')
        #     print('                                                                                     ')
        #     print('                                                                                     ')
        # else:
        #     exit_cleanup_fn = stage.cleanup_fn
        #     input_dir = os.path.join(dirs['stage_dir'], 'input')
        #     scripts_dir = os.path.join(dirs['stage_dir'], 'scripts')
        #     exit_cleanup_fn(run_dir, dirs['stage_dir'], dirs['temp_dir'], scripts_dir, input_dir, batch_size)    
        print('                                                                                     ')
        print('                                                                                     ')
        print(f'All jobs have successfully completed. Exiting {stage.name}.                         ')
        print('                                                                                     ')
        print('                                                                                     ')
        print('=====================================================================================')
        print('=====================================================================================')

    print('                                                                                     ')
    print('                                                                                     ')
    print('                                                                                     ')
    print('                                                                                     ')
    print('=====================================================================================')
    print('=====================================================================================')
    print('                                                                                     ')
    # shutil.copy('/cbica/home/grahamth/sweetnothings/sweetnothings_opencl.py', run_dir)
    # home_base_dir = os.path.join(PIPELINE.home_base, PIPELINE.name)
    # os.mkdir(home_base_dir)
    # copytree(run_dir, home_base_dir, copy_function=copy2_verbose)
    print('                                                                                     ')
    print(f'Run complete. All stages have successfully exited.                                                                      ')
    print('=====================================================================================')
    print('=====================================================================================')




        


    #print(list_files(config.WORKING_DIR))

if __name__ == '__main__':
    run()
