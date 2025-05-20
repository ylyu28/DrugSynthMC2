### Receptor preparation
# receptor pdb files is expected to be extracted by the user prior to running DrugSynthMC2

import shutil
import subprocess 
import os


prot_ids = [f'clusters_{n}' for n in range(1,6)]

for prot_id in prot_ids:
    source_path = f"/Users/yalilyu/Desktop/MRes_P2_DrugSynthMC2/ar_charmm-gui-4344145477/gromacs/{prot_id}.pdb"
    destination_path = f"/Users/yalilyu/Desktop/code/DrugSynthMC2/docking/ar/ar_structures/AR_{prot_id}.pdb"
    shutil.copyfile(source_path, destination_path)

# for prot_id in prot_ids:
#     source_path = f"/Users/yalilyu/Desktop/MRes_P2_DrugSynthMC2/her2_charmm-gui-4567121597/gromacs/{prot_id}.pdb"
#     destination_path = f"/Users/yalilyu/Desktop/python_test_docking/HER2_clusters/HER2_{prot_id}.pdb"
#     shutil.copyfile(source_path, destination_path)

def pdb_to_pdbqt(input_dir, output_dir, protein, num_files):
    """
    Batch conversion of pdb files to pdbqt files.
    """
    for i in range(1, num_files+1):
        input_file = os.path.join(input_dir, f"{protein}_clusters_{i}.pdb")
        output_file = os.path.join(output_dir, f"{protein}_clusters_{i}.pdbqt")

        command = [
            "obabel",
            input_file,
            "-xr", "-xn", "-xp",
            "-O", output_file
        ]

        subprocess.run(command, check=True)
        print(f"Converted {input_file} to {output_file}")

# Convert AR 5 clusters pdb to pdbqt
pdb_to_pdbqt("/Users/yalilyu/Desktop/code/DrugSynthMC2/docking/ar/ar_structures", "/Users/yalilyu/Desktop/code/DrugSynthMC2/docking/ar/ar_structures", "AR", 5)


