import os, shutil
from pathlib import Path
import pathlib
import collections


# def cleanup_stage_dir(stage_dir, temp_dir, scripts_dir, input_dir, batch_size):
#     cleanup_scripts_dir(stage_dir, scripts_dir, batch_size)
#     cleanup_dir(input_dir)
#     cleanup_dir(temp_dir)


def cleanup_scripts_dir(run_dir, stage_dir, scripts_dir, batch_size):
    print(f'Archiving submission scripts and removing {scripts_dir}')
    script_archive_dir = create_submit_script_archive_dir(run_dir, stage_dir)
    print('Preparing and moving submission scripts...')
    combine_submit_scripts(scripts_dir, script_archive_dir)
    if collections.Counter(p.suffix for p in pathlib.Path(script_archive_dir).iterdir())['.sh'] == batch_size:
        scripts_archive = make_directory_archive(script_archive_dir)
        print(f'Submission scripts were archived to: {scripts_archive}')
        remove_dir_tree(script_archive_dir)
        remove_dir_tree(scripts_dir)
    else:
        print(f'WARNING there is a mismatch between number of batches and submit scripts, the scripts dir was not removed.')
        scripts_archive = make_directory_archive(script_archive_dir)
        print(f'Submission scripts were archived to: {scripts_archive}')

def create_submit_script_archive_dir(run_dir, stage_dir):
    script_archive_dir = os.path.join(run_dir, pathlib.Path(stage_dir).stem+'_submit_script_archive')
    os.mkdir(script_archive_dir)
    return script_archive_dir


def combine_submit_scripts(scripts_dir, script_archive_dir):
    rootdir = scripts_dir
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            #print os.path.join(subdir, file)
            filepath = Path(subdir + os.sep + file)
#             print(filepath)
            if filepath.stem == 'submit_jobs':
                parent_dir = filepath.parent
                batched_script_path = os.path.join(script_archive_dir, '%s.%s.sh'%(parent_dir.stem,filepath.name))
                shutil.move(filepath, batched_script_path)  

def cleanup_input_dir(run_dir, stage_dir, input_dir, batch_size):
    print(f'Archiving sgeout logs and removing {input_dir}')
    sgeout_archive_dir = create_sgeout_archive_dir(run_dir, stage_dir)
    print('Preparing and moving sgeout logs...')
    combine_sgeout_logs(input_dir, sgeout_archive_dir)
    sgeout_archive = make_directory_archive(sgeout_archive_dir)
    print(f'Submission scripts were archived to: {sgeout_archive}')
    remove_dir_tree(sgeout_archive_dir)
    remove_dir_tree(input_dir)



def create_sgeout_archive_dir(run_dir, stage_dir):
    sgeout_archive_dir = os.path.join(run_dir, pathlib.Path(stage_dir).stem+'_sgeout_archive')
    os.mkdir(sgeout_archive_dir)
    return sgeout_archive_dir

def combine_sgeout_logs(scripts_dir, sgeout_archive_dir):
    rootdir = scripts_dir
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            #print os.path.join(subdir, file)
            filepath = Path(subdir + os.sep + file)
#             print(filepath)
            if filepath.stem == 'submit_jobs.sh':
                parent_dir = filepath.parent
                batched_script_path = os.path.join(sgeout_archive_dir, '%s.%s.log'%(parent_dir.stem,filepath.name))
                shutil.move(filepath, batched_script_path) 


def cleanup_dir(dir_path):
    remove_dir_tree(dir_path)
    print(f'{dir_path} was removed.')

def remove_dir_tree(target_dir):
    try:
        shutil.rmtree(target_dir)
    except OSError as e:
        print("Error: %s : %s" % (target_dir, e.strerror))

def make_directory_archive(source):
    destination = os.path.join(source + '.tar.gz') 
    # format = 'gztar'

    base = os.path.basename(destination)
    # print(base)
    name = base.split('.')[0]
    # print(name)
    format = 'gztar'
    zip_ext = 'tar.gz'
    # print(format)
    archive_from = os.path.dirname(source)
    # print(archive_from)
    archive_to = os.path.basename(source.strip(os.sep))
    # print(archive_to)

    shutil.make_archive(name, format, archive_from, archive_to)
    shutil.move('%s.%s'%(name,zip_ext), destination)
    return destination

