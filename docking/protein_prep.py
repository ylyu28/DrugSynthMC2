### Receptor preparation
# receptor pdb files is expected to be extracted by the user prior to running DrugSynthMC2

import shutil
import subprocess 
import os


# def pdb_to_pdbqt(input_dir, output_dir, protein, num_files):
#     """
#     Batch conversion of pdb files to pdbqt files.
#     """
#     for i in range(1, num_files+1):
#         input_file = os.path.join(input_dir, f"{protein}_clusters_{i}.pdb")
#         output_file = os.path.join(output_dir, f"{protein}_clusters_{i}.pdbqt")

#         command = [
#             "obabel",
#             input_file,
#             "-xr", "-p 7.4", "--partialcharge eem",
#             "-O", output_file
#         ]

#         subprocess.run(command, check=True)
#         print(f"Converted {input_file} to {output_file}")

# # Convert AR 4 clusters pdb to pdbqt
# pdb_to_pdbqt("./docking/ar/ar_structures", "./docking/ar/ar_structures", "AR", 4)

# Used cluster obabel to convert pdb to pdbqt
# Command line:
# obabel ar_clusters_1.pdb -O ar_clusters_1.pdbqt -xr -p 7.4 --partialcharge eem 
