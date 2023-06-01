import os
from tkinter import filedialog
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
import rmsd
import shutil
import warnings
from Bio import BiopythonWarning
import pymol
import time


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
        atoms = []
        for residue in model.get_residues():
            res_id = residue.get_id()[1]
            if res_id not in missing_residues_all:
                if (min_res <= res_id) and (res_id <= max_res):
                    for atom_name in ['CA', 'CB']:
                        try:
                            atom = residue[atom_name]
                            atoms.append(atom)
                        except KeyError:
                            pass
                    # ca_atoms.append(residue["CA"])
                    # except KeyError:
                    #     pass
        coords.append(np.array([atom.get_coord() for atom in atoms]))
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

def get_rmsd_matrix(coords):
    n = len(coords)
    rmsd_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            rmsd_matrix[i][j] = rmsd_matrix[j][i] = rmsd.kabsch_rmsd(coords[i], coords[j])
    return rmsd_matrix

def hierarchical_clustering(rmsd_matrix, method='average', filename='dendrogram.png', labels=None):
    linkage_matrix = linkage(pdist(rmsd_matrix), method=method)
    fig, ax = plt.subplots(figsize=(24, 15))
    dn = dendrogram(linkage_matrix, ax=ax, labels=labels, orientation='left', color_threshold=1)
    categorized_label = {}
    for i, label_idx in enumerate(dn['leaves']):
        label = labels[label_idx]
        color = dn['leaves_color_list'][i]
        categorized_label[label] = color

    print(categorized_label)
    plt.savefig(os.path.join(dir_path,filename))
    return categorized_label

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

def Hcluster(dir_path):
    dir_path = dir_path
    files = [os.path.join(dir_path, f) for f in os.listdir(dir_path) if f.endswith(".pdb")]

    labels = [os.path.splitext(os.path.basename(file))[0] for file in files]
    # labels = [label.split('_')[0] + '_chain_' + label.split('_')[2] for label in labels]

    min_res, max_res = get_min_max_residues(files)
    missing_residues_all = remove_missing_residues_all(files, min_res, max_res)
    coords = get_coords(files, min_res, max_res, missing_residues_all, dir_path)
    rmsd_matrix = get_rmsd_matrix(coords)
    categorized_label = hierarchical_clustering(rmsd_matrix, labels=labels)

    # Start PyMOL and load the reference pdb file
    pymol.finish_launching()
    time.sleep(1)  # wait for PyMOL to fully load

    color_dict = {'C0': 'blue', 'C1': 'orange', 'C2': 'green', 'C3': 'red', 'C4': 'purple', 'C5': 'brown', 'C6': 'pink'
                  , 'C7': 'gray', 'C8': 'yellow', 'C9': 'cyan'}

    for pdb_file in os.listdir(dir_path):
        if pdb_file.endswith(".pdb"):

            # Load the pdb file
            pymol.cmd.load(os.path.join(dir_path, pdb_file), pdb_file[:-4])
            # time.sleep()
            color_code = categorized_label[pdb_file[:-4]]
            color = color_dict[color_code]
            pymol.cmd.color(color, f"{pdb_file[:-4]}")
            pymol.cmd.set('cartoon_fancy_helices', 1)
            pymol.cmd.set('cartoon_oval_length', 1.1)

    new_session_name = f"aligned.pse"
    pymol.cmd.save(os.path.join(dir_path, new_session_name))
    pymol.cmd.sync()  # wait for the save command to complete
    pymol.cmd.delete('all')

    # Quit PyMOL
    pymol.cmd.quit()

dir_path = filedialog.askdirectory(title='Select a directory to work with')
Hcluster(dir_path)


