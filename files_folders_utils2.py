import os
import shutil
import glob
import config
import time
import json
from functools import partial
import config
import subprocess
import csv
from ast import literal_eval
import shlex
import re
import collections
from rdkit import Chem


def compound_indexer(target_dir, subjob_input_db_ext):
    csv_columns = ['index_id', 'compound_name']
    with open(f'{target_dir}/compound_index_key.csv', 'w') as csvfile:
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
                    d = dict();
                    d['index_id'] = index_id
                    d['compound_name'] = compound_name
                    writer.writerow(d)



def move_files_between_stages(stage_count, run_dir, current_stage_temp_dir, subjob_input_db_ext):
    previous_stage_dir = os.path.join(run_dir, f'stage_{stage_count}')
    target_dir = current_stage_temp_dir
    source_dir = os.path.join(previous_stage_dir, 'output')
    l_of_directories = glob.glob(os.path.join(source_dir, "*", ""))
    # print(f'The list of directories in previous output are: {l_of_directories}')
    for d in range(len(l_of_directories)):
        batch_dir = l_of_directories[d]
        # print(f'The batch_dir is: {batch_dir}')
        file_names = glob.glob(os.path.join(batch_dir, "*" + subjob_input_db_ext))
        # print(f'The file names are: {file_names}')
        for file_name in file_names:
            shutil.copy(os.path.join(file_name), target_dir)

def make_temp_dir(stage_dir):
    temp_dir_path = os.path.join(stage_dir, 'temp')
    temp_dir = os.mkdir(temp_dir_path)
    return temp_dir_path


def setup_run_dir(working_dir, job_name):
    job_path = os.path.join(working_dir, job_name)
    if os.path.exists(job_path):
        print('WARNING: Job directory already exists - exiting.')
        return None
    os.mkdir(job_path)
    return job_path

def make_stage_dir(run_dir, stage):
    job_path = os.path.join(run_dir, stage)
    if os.path.exists(job_path):
        print('WARNING: Job directory already exists - exiting.')
        return None
    os.mkdir(job_path)
    return job_path

def makefolders(root_dir, subfolders):
    concat_path = partial(os.path.join, root_dir)
    makedirs = partial(os.makedirs, exist_ok=True)  # Python 3.2+
    for path in map(concat_path, subfolders):
        makedirs(path)


def list_files(startpath):
    for root, dirs, files in os.walk(startpath):
        level = root.replace(startpath, '').count(os.sep)
        indent = ' ' * 4 * (level)
        print('{}{}/'.format(indent, os.path.basename(root)))
        subindent = ' ' * 4 * (level + 1)
        for f in files:
            print('{}{}'.format(subindent, f))

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def script_maker(input_prefix, batch_number, stage_dir, input_db_ext, output_db_ext, program_exe):
    
    input_dir = os.path.join(stage_dir, 'input', batch_number)
    scripts_dir = os.path.join(stage_dir, 'scripts', batch_number)
    output_dir = os.path.join(stage_dir, 'output', batch_number)
    
    input_file_path = os.path.join(input_dir, input_prefix + input_db_ext)
    
    job_filename = f'run_{input_prefix}.sh'
    job_file_path = os.path.join(scripts_dir, job_filename)
    
    output_filename = f'{input_prefix}{output_db_ext}'
    output_file_path = os.path.join(output_dir, output_filename)

    #print(f'the input file path is {input_file_path}')
    #print(f'the output file path is {output_file_path}')
    #print(f'the jobs file path is {job_file_path}')

    with open(job_file_path, 'w') as job_f:

        exe_cmd = program_exe.format(input_file_path=input_file_path, output_file_path=output_file_path)
        shell_script = '''\
            {exe_cmd}
        '''.format(exe_cmd=exe_cmd)
        job_f.write(shell_script)

def master_chunk_depositor(dir_of_folders, dir_of_files, files_per_folder, input_db_ext, stage_dir, exe_cmd, output_db_ext): 
    
    dir1 = dir_of_files
    list_of_files = glob.glob(os.path.join(dir1, "*" + input_db_ext))
    json_list_of_files = os.path.join(dir1 + 'file_index.json')
    with open(json_list_of_files, 'w') as fp:
        json.dump(list_of_files, fp)

    dir2 = dir_of_folders
    list_of_directories = glob.glob(os.path.join(dir2, "*", ""))
    json_list_of_folders = os.path.join(dir2 + 'directory_index.json')
    with open(json_list_of_folders, 'w') as fp:
        json.dump(list_of_directories, fp)
        
    with open(json_list_of_folders) as folders_in:
        lis=json.load(folders_in)
    with open(json_list_of_files) as files_in:
        lisb=json.load(files_in)
    
    n = files_per_folder
    
    chunked_lisb = chunks(lisb, n)
    lisb_f = dict(enumerate(chunked_lisb))

    i = 0
    j = 0
    m = 0
    a = 0

    for m in range(len(lis)):
        #print(f'total number in lis (m) {len(lis)}')
        #print(f'The variable m is = {m}')
        for j in range(len(lisb_f)):
            #print(f'total number in lisb_f (j) {len(lisb_f)}') 
            sublist = lisb_f[j]
            for a in range(len(sublist)):
                    #print(f'total number in sublist (a) {len(sublist)}')
                    #print(sublist[i+a], lis[m])
                    shutil.copy(sublist[i+a], lis[m])
                    input_f = (sublist[i+a])
                    input_prefix = input_f.split('/')[-1].split('.')[0]
                    # print(f'the input prefix is: {input_prefix}')
                    batch_path = lis[m]
                    batch_number = os.path.basename(os.path.normpath(batch_path))
                    # print(f'the batch number is: {batch_number}')
                    script_maker(input_prefix, batch_number, stage_dir, input_db_ext, output_db_ext, exe_cmd)
            m += 1
        if m >= len(lis):
            break



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



def make_submit_script(stage_dir, scripts_batch_directories):
    
    for batch_directories in scripts_batch_directories:
        batch_directory_path = os.path.join(stage_dir, batch_directories)
        job_files = os.listdir(batch_directory_path)
        file_path = os.path.join(batch_directory_path, 'submit_jobs.sh')
        print(batch_directory_path)
        with open(file_path, 'w') as submit_f:
            for job_f in job_files:
                job_file_path = os.path.join(batch_directory_path, job_f)
                submit_cmd = f'qsub {job_file_path} &\n'
                submit_f.write(submit_cmd)
        print(file_path)

def make_batch_submit_script(stage_dir, scripts_batch_directories, cpu_or_gpu, job_dependency, input_db_ext, run_dir):
    csv_columns = ['run', 'stage', 'batch', 'inputbatchpath', 'numberoffilesin', 'sgecommand', 'sgejobnumber']
    # # csv_file = print(f'{print(os.path.basename(run_dir))}.csv')
    csv_file = os.path.join(run_dir, os.path.basename(stage_dir) + '_submission_metadata.csv')
    with open(csv_file, 'w') as csvfile:
        
        writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
        writer.writeheader()
        for batch_directories in scripts_batch_directories:
            batch_directory_path = os.path.join(stage_dir, batch_directories)
            read_files = glob.glob(os.path.join(batch_directory_path, "*.sh"))
            file_path = os.path.join(batch_directory_path, 'submit_jobs.sh')
            #print(batch_directory_path)
            with open(file_path, "wb") as outfile:
                for f in read_files:
                    with open(f, "rb") as infile:
                        outfile.write(infile.read())
            #print(file_path)
            
            ee = add_job_parameters(file_path, batch_directory_path, cpu_or_gpu, stage_dir, job_dependency, input_db_ext, run_dir)
        
            writer.writerow(ee)



def submit_sge_scripts(cpu_or_gpu, input_batch_path, sge_ready_script_path, input_db_ext, stage_dir, run_dir):
    submit_cmd = [f'{cpu_or_gpu} -wd {input_batch_path} {sge_ready_script_path}']
    
    t = subprocess.check_output(submit_cmd, shell=True)
    x = [int(s) for s in t.split() if s.isdigit()]
    job_number = str(x)[1:-1] 



    batch_number = os.path.basename(input_batch_path)
    job_subprocess_count = len(glob.glob1(input_batch_path, (f'*{input_db_ext}')))
    run_name = os.path.basename(run_dir)
    stage_number = os.path.basename(stage_dir)
    d = dict();
    d['run'] = run_name
    d['stage'] = stage_number
    d['batch'] = batch_number
    d['inputbatchpath'] = input_batch_path
    d['numberoffilesin'] = job_subprocess_count
    d['sgecommand'] = cpu_or_gpu
    d['sgejobnumber'] = job_number
    # print(d)

    return d


def add_job_parameters(sge_ready_script_path, batch_directory_path, cpu_or_gpu, stage_dir, job_dependency, input_db_ext, run_dir):
        # print(sge_ready_script_path)
        f = open(sge_ready_script_path,'r+')
        lines = f.readlines() # read old content
        f.seek(0) # go back to the beginning of the file
        f.write(job_dependency) # write new content at the beginning
        for line in lines: # write old content after new
            f.write(line)
        f.close()
        batch_number = os.path.basename(os.path.normpath(batch_directory_path))
        input_batch_path = os.path.join(stage_dir, 'input', batch_number)
        job_submit_info = submit_sge_scripts(cpu_or_gpu, input_batch_path, sge_ready_script_path, input_db_ext, stage_dir, run_dir)
        # os.system(f'{cpu_or_gpu} -wd {input_batch_path} {sge_ready_script_path}')
        # print(f'The submit command to be executed is: {cpu_or_gpu} -wd {input_batch_path} {sge_ready_script_path}')
        # print(f'{sge_ready_script_path} submitted as an SGE job submitted successfully.')
        return job_submit_info

def check_job_status(run_dir, stage_dir):
    csv_file = os.path.join(run_dir, os.path.basename(stage_dir) + '_submission_metadata.csv')
    # open the file in universal line ending mode 
    with open(csv_file, 'rU') as infile:
      # read the file as a dictionary for each row ({header : value})
      reader = csv.DictReader(infile)
      data = {}
      for row in reader:
        for header, value in row.items():
          try:
            data[header].append(value)
          except KeyError:
            data[header] = [value]

    # extract the variables you want
    names = data['sgejobnumber']
    return names

def trackJobs(jobs, stage_dir, batch_size, input_batch_directories, output_batch_directories, input_db_ext, output_db_ext, waittime=15):
    
    while len(jobs) != 0:
        for jobid in jobs:
            # print(len(jobs))
            # print(jobid)
            x = subprocess.Popen(['qstat', '-j', jobid], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            std_out, std_err = x.communicate()
            if std_err :
                jobs.remove(jobid)
                print(f'{jobid} is now complete.')
                # break
        input_vs_output_count(stage_dir, batch_size, input_batch_directories, output_batch_directories, input_db_ext, output_db_ext)
        os.system("sleep " + str(waittime))
        print(f'There are {len(jobs)} jobs remaing. Waiting {str(waittime)}s until checking again for job completion.')
        print(f'The jobs remaining are {(jobs)}')

def file_counter(list_of_directories, file_ext):
    file_count = []
    for batch_directory in list_of_directories:
        # print(list_of_directories)
        number_of_files = len([f for f in os.listdir(batch_directory) if re.search(rf'({file_ext})$', f)])
        file_count.append(number_of_files)

    return file_count

def sum_files_in_count(file_count):
    num_file_count = file_count
    for idx, file in enumerate(num_file_count):
        num_file_count[idx] = int(file)
    
    total_count = sum(num_file_count)

    return total_count

def input_vs_output_count(stage_dir, batch_size, input_batch_directories, output_batch_directories, input_db_ext, output_db_ext): 

    # input_batch_directories = ([r'input/batch{:04d}'.format(i) for i in range(1, batch_size + 1)])
    # output_batch_directories = ([r'output/batch{:04d}'.format(i) for i in range(1, batch_size + 1)])
    number_of_input_files = [os.path.join(stage_dir, x) for x in input_batch_directories]
    number_of_output_files = [os.path.join(stage_dir, x) for x in output_batch_directories]
    # print(number_of_input_files)
    # print(number_of_output_files)
    # while complete_perc <= 100:
    if file_counter(number_of_input_files, input_db_ext) is not None:
        file_count_input = file_counter(number_of_input_files, input_db_ext)
        total_input = sum(file_count_input)
        # print(file_counter(number_of_input_files, input_db_ext))
        print(f'The total number of input files is: {total_input}')

    if file_counter(number_of_output_files, output_db_ext) is not None:
        file_count_output = file_counter(number_of_output_files, output_db_ext)
        print(file_count_output)
        total_output = sum(file_count_output)
        # print(file_counter(number_of_output_files, output_db_ext))
        print(f'The total number of output files is: {total_output}')
        complete_perc = total_output*100 / total_input
        print(f'The stage is {complete_perc}% complete. All jobs must exit prior to the next stage starting.')

    else:
        print('No output files detected. Checking again momentarily.')

def resubmit_failed_jobs(output_batch_directories, stage_dir, run_dir, input_db_ext, output_db_ext):
    csv_columns = ['run', 'stage', 'batch', 'inputbatchpath', 'numberoffilesin', 'sgecommand', 'sgejobnumber']
    csv_file = os.path.join(run_dir, os.path.basename(stage_dir) + '_resubmission_metadata.csv')
    previous_submissioncsv = os.path.join(run_dir, os.path.basename(stage_dir) + '_submission_metadata.csv')
    with open(csv_file, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
        writer.writeheader()
        data = csv_file
        with open(data, 'rU') as infile:
            # read the file as a dictionary for each row ({header : value})
            reader = csv.DictReader(infile)
            data = {}
            for row in reader:
                for header, value in row.items():
                    try:
                        data[header].append(value)
                    except KeyError:
                        data[header] = [value]

            # extract the variables you want

        

        number_of_output_files = [os.path.join(stage_dir, x) for x in output_batch_directories]
        for batch_directory in number_of_output_files:
            number_of_files = len([f for f in os.listdir(batch_directory) if re.search(rf'({output_db_ext})$', f)])
            if number_of_files <= 1:
                with open(previous_submissioncsv, 'rU') as infile:
            # read the file as a dictionary for each row ({header : value})
                    reader = csv.DictReader(infile)
                    dataold = {}
                    for row in reader:
                        for header, value in row.items():
                            try:
                                dataold[header].append(value)
                            except KeyError:
                                dataold[header] = [value]
                print(batch_directory)
                print(f'The following batch has a failed job: {batch_directory}')
                batch_num = os.path.basename(batch_directory)
                print(batch_num)
                isolate_batch_num = re.findall("\d+", batch_num)[0]
                        
                       
                clean = int(isolate_batch_num)
                key = clean - 1
                
                print(clean)
                print(key)
                input_path = os.path.join(stage_dir, 'input', batch_directory)
                script_path = os.path.join(stage_dir, 'scripts', batch_num, 'submit_jobs.sh')
                print(script_path)
                sge_ready_script_path = script_path
                cpu_or_gpu = dataold['sgecommand'][key]
                print(cpu_or_gpu)
                input_batch_path = input_path
                
                failed_entry = submit_sge_scripts(cpu_or_gpu, input_batch_path, sge_ready_script_path, input_db_ext, stage_dir, run_dir)
                writer.writerow(failed_entry)
            else:
                print(f'No failed job was detected in: {batch_directory}')

def check_failed_job_status(run_dir, stage_dir):
    csv_file = os.path.join(run_dir, os.path.basename(stage_dir) + '_resubmission_metadata.csv')
    # open the file in universal line ending mode 
    with open(csv_file, 'rU') as infile:
      # read the file as a dictionary for each row ({header : value})
      reader = csv.DictReader(infile)
      data = {}
      for row in reader:
        for header, value in row.items():
          try:
            data[header].append(value)
          except KeyError:
            data[header] = [value]

    # extract the variables you want
    failed_jobs_nums = data.get('sgejobnumber')
    if failed_jobs_nums:
        failed_jobs_list = failed_jobs_nums
    else:
        failed_jobs_list = 0
    print(failed_jobs_list)
    return failed_jobs_list

# def check_sge_job_list
#         for jobid in jobid_list_of_submitted_jobs
#             submit_cmd = 'qstat'
#             str = subprocess.check_output(qstat_cmd, shell=True)
#             match_pid = blah blah
#             if match_pid is True
#                 remove from list
#             else
#                 continue


# def add_job_parameters_nodependency(sge_ready_script_path, batch_directory_path, cpu_or_gpu):
#         print(sge_ready_script_path)
#         # f = open(sge_ready_script_path,'r+')
#         # lines = f.readlines() # read old content
#         # f.seek(0) # go back to the beginning of the file
#         # f.write(nvidia-smi) # write new content at the beginning
#         # for line in lines: # write old content after new
#         #     f.write(line)
#         # f.close()
#         batch_number = os.path.basename(os.path.normpath(batch_directory_path))
#         input_batch_path = os.path.join(stage_dir, 'input', batch_number)
#         os.system(f'{cpu_or_gpu} -wd {input_batch_path} {sge_ready_script_path}')
#         os.system(f'{cpu_or_gpu} -wd {batch_directory_path} {sge_ready_script_path}')
#         print(f'The submit command to be executed is: {cpu_or_gpu} -wd {batch_directory_path} {sge_ready_script_path}')
#         print(f'{sge_ready_script_path} submitted as an SGE job submitted successfully.')
#         time.sleep(0.2)

    #scripts_dir = os.path.join(job_path, 'scripts')
    #input_dir = os.path.join(job_path, 'input')
    #output_dir = os.path.join(job_path, 'output')
    
    #os.mkdir(scripts_dir)
    #os.mkdir(input_dir)
    #os.mkdir(output_dir)
    #return job_path

# from __future__ import print_function
# import os
# import sys

# def find_dirs(root_dir, names, exclude_folders=[]):
#     try:
#         for entry in os.listdir(root_dir):
#             entry_path = os.path.join(root_dir, entry)
#             entry_lowercase = entry.lower()
#             if os.path.isdir(entry_path):
#                 if entry_lowercase in names:
#                     yield entry_path
#                 elif entry_lowercase not in exclude_folders:
#                     for result in find_dirs(entry_path, names, exclude_folders):
#                         yield result
#     except OSError as e:
#         print(e, file=sys.stderr)

# product_dirs = []
# event_log_path = '/tmp/findlog.txt'
# with open(event_log_path, 'a') as log:
#     for lib in find_dirs('/', ['lib', 'library'], ['user']):
#         print(lib)
#         product_dirs.append(lib)
#         log.write(lib + "\n")
#         break      # Stop after finding the first match