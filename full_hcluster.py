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
from consv_to_pdb import readfile

def get_resi_atoms(files):
    common_resi_atoms = []
    for file in files:
        with open(file, 'r') as f:
            lines = f.readlines()

        resi_atoms = []

        for line in lines:
            if line.startswith('ATOM'):
                resi_atom = line[22:26].strip(' ') + '_' + line[13:16].strip(' ')
                resi_atoms.append(resi_atom)
        if not common_resi_atoms:
            common_resi_atoms = resi_atoms
        common_resi_atoms = list(set(common_resi_atoms).intersection(resi_atoms))

    common_resi_atoms = sorted(common_resi_atoms)
    print(common_resi_atoms)
    resi_number = list(set([i.split('_')[0] for i in common_resi_atoms]))
    print(len(common_resi_atoms))
    print(resi_number)
    print(len(resi_number))
    return common_resi_atoms

def get_resi_atoms_all(files, seqs_dictionary):
    common_resi_atoms = []
    for file in files:
        accession = os.path.splitext(os.path.basename(file))[0].split('_')[1]
        with open(file, 'r') as f:
            lines = f.readlines()

        resi_atoms = []

        for line in lines:
            if line.startswith('ATOM') and int(line[22:26].strip(' ')) in seqs_dictionary[accession][1].keys():
                resi_atom = str(seqs_dictionary[accession][1][int(line[22:26].strip(' '))]) + '_' + line[13:16].strip(' ')
                resi_atoms.append(resi_atom)
        if not common_resi_atoms:
            common_resi_atoms = resi_atoms
        common_resi_atoms = list(set(common_resi_atoms).intersection(resi_atoms))

    common_resi_atoms = sorted(common_resi_atoms)
    resi_number = sorted(list(set([i.split('_')[0] for i in common_resi_atoms])))
    print(common_resi_atoms)
    print(len(common_resi_atoms))
    print(resi_number)
    print(len(resi_number))
    return common_resi_atoms


# def get_coords(files, common_resi_atoms, dir_path):
#     warnings.simplefilter('ignore', BiopythonWarning)
#     coords = []
#     parser = PDBParser()
#     for file in files:
#         # Copy the original file to a temporary file
#         temp_file = os.path.join(dir_path, f"temp_{os.path.basename(file)}")
#         shutil.copy2(file, temp_file)
#         structure = parser.get_structure(os.path.basename(file), temp_file)
#         model = structure[0]
#         atoms = []
#
#         for residue in model.get_residues():
#             res_id = residue.get_id()[1]
#             for resi_atoms in common_resi_atoms:
#                 if int(resi_atoms.split('_')[0]) == res_id:
#                     atoms.append(residue[resi_atoms.split('_')[1]])
#         coords.append(np.array([atom.get_coord() for atom in atoms]))
#         # Delete the temporary file
#         os.remove(temp_file)
#
#     print(np.array(coords).shape)
#     return coords

def get_coords(files, common_resi_atoms, dir_path):
    warnings.simplefilter('ignore', BiopythonWarning)
    coords = []
    parser = PDBParser()

    # Create a 'temp' folder in dir_path if it doesn't exist
    temp_folder = os.path.join(dir_path, 'temp')
    os.makedirs(temp_folder, exist_ok=True)

    for file in files:
        # Copy the original file to a temporary file in the 'temp' folder
        temp_file = os.path.join(temp_folder, f"temp_{os.path.basename(file)}")
        shutil.copy2(file, temp_file)

        structure = parser.get_structure(os.path.basename(file), temp_file)
        model = structure[0]
        atoms = []

        for residue in model.get_residues():
            res_id = residue.get_id()[1]
            for resi_atoms in common_resi_atoms:
                if int(resi_atoms.split('_')[0]) == res_id:
                    atoms.append(residue[resi_atoms.split('_')[1]])

        coords.append(np.array([atom.get_coord() for atom in atoms]))

        # Delete the temporary file
        os.remove(temp_file)

    # Delete the 'temp' folder
    os.rmdir(temp_folder)

    return coords



def get_coords_all(files, common_resi_atoms, dir_path, seqs_dictionary):
    warnings.simplefilter('ignore', BiopythonWarning)
    coords = []
    parser = PDBParser()
    for file in files:
        accession = os.path.splitext(os.path.basename(file))[0].split('_')[1]
        # Copy the original file to a temporary file
        temp_file = os.path.join(dir_path, f"temp_{os.path.basename(file)}")
        shutil.copy2(file, temp_file)
        structure = parser.get_structure(os.path.basename(file), temp_file)
        model = structure[0]
        atoms = []

        for residue in model.get_residues():
            res_id = residue.get_id()[1]
            for resi_atoms in common_resi_atoms:
                if seqs_dictionary[accession][2][int(resi_atoms.split('_')[0])] == res_id:
                    atoms.append(residue[resi_atoms.split('_')[1]])
        coords.append(np.array([atom.get_coord() for atom in atoms]))
        # Delete the temporary file
        os.remove(temp_file)

    print(np.array(coords).shape)
    return coords

def get_rmsd_matrix(coords):
    n = len(coords)
    rmsd_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            rmsd_matrix[i][j] = rmsd_matrix[j][i] = rmsd.kabsch_rmsd(coords[i], coords[j])
    return rmsd_matrix

def hierarchical_clustering(rmsd_matrix, method='complete', filename='dendrogram.png', labels=None, dir_path=None):
    linkage_matrix = linkage(pdist(rmsd_matrix), method=method)
    fig, ax = plt.subplots(figsize=(24, 15))
    dn = dendrogram(linkage_matrix, ax=ax, labels=labels, orientation='left', color_threshold=2.5)
    categorized_label = {}
    for i, label_idx in enumerate(dn['leaves']):
        label = labels[label_idx]
        color = dn['leaves_color_list'][i]
        categorized_label[label] = color

    print(categorized_label)
    plt.savefig(os.path.join(dir_path,filename))
    return categorized_label


def H_cluster():
    dir_path = filedialog.askdirectory(title='Select a directory to work with')
    # # ref_pdb_file = filedialog.askopenfilename(title='Select the file you want to use as the reference',initialdir=dir_path)
    files = [os.path.join(dir_path, f) for f in os.listdir(dir_path) if f.endswith(".pdb")]
    labels = [os.path.splitext(os.path.basename(file))[0] for file in files]
    common_resi_atoms = get_resi_atoms(files)
    coords = get_coords(files, common_resi_atoms, dir_path)
    rmsd_matrix = get_rmsd_matrix(coords)
    categorized_label = hierarchical_clustering(rmsd_matrix, labels=labels, dir_path=dir_path)

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

def H_cluster_all_type():


    dir_path = filedialog.askdirectory(title='Select a directory to work with')
    filename = filedialog.askopenfilename(title='Select the fasta file')
    seqs, headers, seqs_dictionary, missing_map = readfile(filename)
    files = [os.path.join(dir_path, f) for f in os.listdir(dir_path) if f.endswith(".pdb")]
    labels = [os.path.splitext(os.path.basename(file))[0] for file in files]
    print(files)
    print(labels)
    common_resi_atoms = get_resi_atoms_all(files, seqs_dictionary)
    coords = get_coords_all(files, common_resi_atoms, dir_path, seqs_dictionary)
    rmsd_matrix = get_rmsd_matrix(coords)
    categorized_label = hierarchical_clustering(rmsd_matrix, labels=labels, dir_path=dir_path)

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






########################################
# H_cluster()
H_cluster_all_type()