import os
import time
from tkinter import filedialog
from tkinter import Tk
import pymol


root = Tk()
root.withdraw()

def align_and_save_all_pdb_files(directory, ref_pdb_file):
    # Create a new directory for aligned pdb files
    parent_dir = os.path.dirname(directory)
    new_dir_name = os.path.join(parent_dir, f"{os.path.basename(directory)}_aligned")
    os.makedirs(new_dir_name, exist_ok=True)

    # Start PyMOL and load the reference pdb file
    pymol.finish_launching()
    time.sleep(1)  # wait for PyMOL to fully load
    pymol.cmd.load(ref_pdb_file, 'ref')
    distance_cutoff = 3.0
    save = 1

    # Loop through all pdb files in the directory
    for pdb_file in os.listdir(directory):
        if pdb_file.endswith(".pdb"):

            # Load the pdb file
            pymol.cmd.load(os.path.join(directory, pdb_file), pdb_file[:-4])

            # Align the pdb file to the reference pdb file
            pymol.cmd.align(pdb_file[:-4], 'ref')

            pymol.cmd.origin('ref')
            pymol.cmd.zoom("all", 0.8)
            #
            # Save the aligned pdb file with a new name
            new_file_name = f"{'_'.join(pdb_file[:-4].split('_')[:-1])}.pdb" #erase_aligned tag
            # new_file_name = f"{'_'.join(pdb_file[:-4].split('_')[:-1])}_aligned.pdb"
            if save != 0:
                pymol.cmd.save(os.path.join(new_dir_name, new_file_name), pdb_file[:-4])
            else: save = 1
            pymol.cmd.sync()  # wait for the save command to complete
            pymol.finish_launching()

            # Delete the loaded pdb file to free up memory
            pymol.cmd.delete(pdb_file[:-4])

    # Delete the reference pdb file to free up memory
    pymol.cmd.delete('ref')

    # Quit PyMOL
    pymol.cmd.quit()

dir_path = filedialog.askdirectory(title='Select a directory to work with')
ref_pdb_file = filedialog.askopenfilename(title='Select the file you want to use as the reference',initialdir=dir_path)
align_and_save_all_pdb_files(dir_path, ref_pdb_file)