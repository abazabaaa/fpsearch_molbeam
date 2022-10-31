import os
import shutil
import glob
from openeye import oechem
import math
import json
import re

def split_molecule_database(input_file, num_parts, stage_dir, temp_dir_path, ext, SDF_name):
    temp_prefix_path = os.path.join(temp_dir_path, SDF_name)
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

def MolCount(input_file):

    moldb = oechem.OEMolDatabase(input_file)
    nummols = moldb.NumMols()
    print("%s contains %d molecule(s)." % (input_file, nummols))
    return nummols

def create_output_stream(outbase, ext, chunk):
    newname = outbase + ('_%05d' % chunk) + ext
    ofs = oechem.oemolostream()
    if not ofs.open(newname):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % newname)
    return ofs

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

# def move_files(temp_dir, stage_dir, input_batch_directories, num_parts, batch_size):
#     ## the total number of files to be distributed into each folder in a list = number of
#     ## molecules (numparts)/()
#     num_files_per_batch = math.floor(int(num_parts)/int(batch_size))

#     files = [os.path.join(temp_dir, f) for f in os.listdir(temp_dir)]

#     list_input_dir = (os.path.join(stage_dir, l) for l in input_batch_directories)

#     #f_base = os.path.basename(f)

#     ##list of lists, chunk of files per batch
#     chunked_files = chunks(files, num_files_per_batch)

#     input_path_listlist = dict(enumerate(chunked_files, 1))
#     listlist_path = os.path.join(temp_dir, 'chunk_list_path_index.json')
#     with open(listlist_path, 'w') as fp:
#     	json.dump(input_path_listlist, fp)

#     with open(listlist_path) as file:
#     	sdf_path_index = json.load(file)
#     	for x in sdf_path_index:
#     		while 

    #with open(index_input_dir) as file2:
    	#input_path_index = json.load(file2)
    	#for y in input_path_index['1']:
    		#print(y)

    #print(list_input_dir)

    # if not os.path.exists(temp_dir):
    #     raise Exception('Directory does not exist ({0}).'.format(temp_dir))

    # i = 0
    # curr_subdir = batch_input_dir

    # for f in files:
    #     # create new subdir if necessary
    #     if i % num_files_per_batch == 0:
    #         subdir_name = os.path.join(stage1/input, 'batch{0:03d}'.format(i / num_files_per_batch + 1))
            
    #         curr_subdir = subdir_name

        # move file to current dir
        
        
        #i += 1