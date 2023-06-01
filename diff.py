import os
from tkinter import filedialog
import numpy as np
from Bio.PDB import PDBParser
import shutil
import warnings
from Bio import BiopythonWarning

def get_coords(files, min_res, max_res, missing_residues_all, dir_path):
    warnings.simplefilter('ignore', BiopythonWarning)
    coords = []
    parser = PDBParser()
    for file in files:
        # Copy the original file to a temporary file
        temp_file = os.path.join(dir_path, f"temp_{os.path.basename(file)}")
        shutil.copy2(file, temp_file)
        structure = parser.get_structure(os.path.basename(file), temp_file)
        model = structure[0]
        ca_atoms = []
        for residue in model.get_residues():
            res_id = residue.get_id()[1]
            if res_id not in missing_residues_all:
                if (min_res <= res_id) and (res_id <= max_res):
                    # try:
                    ca_atoms.append(residue["CA"])
                    # except KeyError:
                    #     pass
        coords.append(np.array([atom.get_coord() for atom in ca_atoms]))
        # Delete the temporary file
        os.remove(temp_file)

    print(np.array(coords).shape)
    return coords


def get_min_max_residues(files):
    first_residues = []
    last_residues = []

    for file in files:
        with open(file, 'r') as f:
            lines = f.readlines()

        first_residue = None
        last_residue = None

        for line in lines:
            if line.startswith('ATOM'):
                residue_id = int(line[22:26].strip())
                if first_residue is None or residue_id < first_residue:
                    first_residue = residue_id
                if last_residue is None or residue_id > last_residue:
                    last_residue = residue_id

        first_residues.append(first_residue)
        last_residues.append(last_residue)

    return max(first_residues), min(last_residues)

def remove_missing_residues(file, min_res, max_res):
    with open(file, 'r') as f:
        lines = f.readlines()

    missing_residues = []

    residue_id = min_res
    for line in lines:
        if line.startswith('ATOM'):
            if (int(line[22:26].strip('')) == residue_id) or (int(line[22:26].strip('')) < residue_id):
                pass
            elif int(line[22:26].strip('')) == residue_id + 1:
                residue_id = residue_id + 1
            else:
                residue_id = residue_id + 1
                missing_residues.append(residue_id)

    return missing_residues

def remove_missing_residues_all(files, min_res, max_res):
    missing_residues_all = []
    for file in files:
        missing_residues = remove_missing_residues(file, min_res, max_res)
        for i in missing_residues:
            missing_residues_all.append(i)

    return list(set(missing_residues_all))


def get_distance(coords, i):
    total_array = coords[0].copy()
    for array in coords[1:]:
        total_array += array
    mean_array = total_array / len(coords)
    differ = coords[i] - mean_array
    seq = [j for j in range(min_res, max_res + 1) if j not in missing_residues_all]
    diff_dict = {}
    total_diff_list = []
    for j in range(len(coords)):
        diff_others = coords[j] - mean_array
        for int in range(len(diff_others)):
            dist = np.sqrt(np.sum(np.square(differ[int])))
            total_diff_list.append(dist)

    # diff_list = [np.sqrt(np.sum(np.square(diff))) for diff in differ]

    for int in range(len(differ)):
        key = seq[int]
        dist = np.sqrt(np.sum(np.square(differ[int])))
        dist_scaled = (dist - min(total_diff_list)) / (max(total_diff_list) - min(total_diff_list))*500
        diff_dict[key] = dist_scaled

    return diff_dict

def get_diffpdb(pdb_file,distance_dict, dir_path):
    diffpdb_dir = os.path.join(dir_path, "diffpdb")
    if not os.path.exists(diffpdb_dir):
        os.makedirs(diffpdb_dir)

    with open(pdb_file, 'r') as f:
        lines = f.readlines()

    for i in range(len(lines)):
        if lines[i].startswith('ATOM'):
            resi_num = int(lines[i][22:26])
            b_factor = float(lines[i][60:66])
            if resi_num in distance_dict.keys():
                distance = distance_dict[resi_num]
            else:
                # If the current resi_num is not in the dictionary, try to calculate the distance
                # based on the previous or next 5 resi_nums
                resi_nums = list(distance_dict.keys())
                prev_resi_list = [x for x in resi_nums if x < resi_num]
                prev_resi_list.sort(reverse=True)
                next_resi_list = [x for x in resi_nums if x > resi_num]
                next_resi_list.sort(reverse=False)
                try:
                    prev_resi_num = prev_resi_list[0]
                except IndexError:
                    prev_resi_num = None
                try:
                    next_resi_num = next_resi_list[0]
                except IndexError:
                    next_resi_num = None
                if prev_resi_num and next_resi_num:
                    # If there are both previous and next resi_nums, interpolate the distance
                    prev_distance = distance_dict[prev_resi_num]
                    next_distance = distance_dict[next_resi_num]
                    distance = prev_distance + (next_distance - prev_distance) * (resi_num - prev_resi_num) / (
                                next_resi_num - prev_resi_num)
                elif prev_resi_num:
                    distance = distance_dict[prev_resi_num] + (distance_dict[prev_resi_num] - distance_dict[prev_resi_list[1]])\
                               *(resi_num - prev_resi_num) / (prev_resi_num - prev_resi_list[1])
                elif next_resi_num:
                    # distance = distance_dict[next_resi_num]
                    distance = distance_dict[next_resi_num] + (distance_dict[next_resi_num] - distance_dict[next_resi_list[1]])\
                               *(next_resi_num - resi_num) / (next_resi_list[1] - next_resi_num)
                else:
                    # If there are no previous or next resi_nums, set the distance to 0
                    distance = 0

            if distance > 500:
                distance = 500
            if distance < 0:
                distance = 0
            lines[i] = lines[i][:60] + f'{distance:6.2f}' + lines[i][66:]

    file_path = os.path.join(diffpdb_dir, f'{os.path.basename(pdb_file)[:-4]}_diff.pdb')
    with open(file_path, 'w') as f:
        f.writelines(lines)

dir_path = filedialog.askdirectory(title='Select a directory to work with')
# ref_pdb_file = filedialog.askopenfilename(title='Select the file you want to use as the reference',initialdir=dir_path)
files = [os.path.join(dir_path, f) for f in os.listdir(dir_path) if f.endswith(".pdb")]


min_res, max_res = get_min_max_residues(files)
missing_residues_all = remove_missing_residues_all(files, min_res, max_res)
print(missing_residues_all)
coords = get_coords(files, min_res, max_res, missing_residues_all, dir_path)
for i in range(len(files)):
    distance_dict = get_distance(coords, i)
    get_diffpdb(files[i], distance_dict, dir_path)






