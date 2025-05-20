try:
    from scrubber import Scrub
except ImportError:
    from molscrub import Scrub
from multiprocessing import Pool
from time import sleep
import os
import math
import subprocess
import time

def prepare_ligand(smile_string):

    """
    molscrub used to assign protonation states and generate ligand conformers
    writer functions in meeko to save the prepared ligand molecules into PDBQT files
    """

    pdbqt_path = "/Users/yalilyu/Desktop/code/DrugSynthMC2/docking/ligand.pdbqt"
    ligand_sdf_path = "/Users/yalilyu/Desktop/code/DrugSynthMC2/docking/ligand.sdf"
    
    with open(pdbqt_path, 'w') as file:
    # erase log contents
        pass
    with open(ligand_sdf_path, 'w') as file:
    # erase log ligand_sdf_path
        pass
    subprocess.run([
        "scrub.py",
        smile_string,
        "-o", ligand_sdf_path,
        "--skip_tautomers",
        "--ph_low", "7",
        "--ph_high", "7"
    ], check=True)

    subprocess.run([
        "mk_prepare_ligand.py",
        "-i", ligand_sdf_path,
        "-o", pdbqt_path
    ], check=True)

    with open(pdbqt_path, 'r') as f:
        ligand_content = f.read()

    return pdbqt_path, ligand_content

def extract_affinity_table_from_stdout(stdout):

    """
    Extract only the affinity table block from the Vina output stdout.
    """

    lines = stdout.splitlines()
    table_started = False
    table_lines = []

    for line in lines:
        if line.strip().startswith("-----+"):
            table_started = True
            table_lines.append(line)
        elif table_started:
            # The table ends after the separator line with dashes
            if line.strip().startswith("-----+"):
                table_lines.append(line)
            elif line.strip() == '' or line.startswith('#'):
                continue
            elif line.startswith('---'):
                # End of table
                break
            else:
                if line.strip() != '':
                    table_lines.append(line)
    return '\n'.join(table_lines)

def run_vina(args):

    print(f"Process {os.getpid()} starting docking for cluster {args[-1]}")

    protein, ligand_pdbqt_path, config_file, protein_file, cluster_index = args
    command = [
    "vina",
    "--receptor", protein_file,
    "--ligand", ligand_pdbqt_path,
    "--config", config_file
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    stdout = result.stdout
    affinities_table = extract_affinity_table_from_stdout(stdout)

    log = f"\n--- Affinity Table for {protein} Cluster {cluster_index} ---\n"
    log += affinities_table
    log += '\n\n'
    
    print(f"Process {os.getpid()} finished docking for cluster {args[-1]}")
    return log



def vina_multiprocessing(protein, ligand_pdbqt_path, config_file, num_files):

    """
    Parallel Docking of ligand to cluster structures of a specified protein using multiprocessing
    and write affinities table to log, then extract the most favorable affinity value.
    """
    
    start_time = time.time()
    log_path = f'/Users/yalilyu/Desktop/code/DrugSynthMC2/docking/{protein}.txt'

    with open(log_path, 'w') as f:
        pass  # create an empty file
    
    args_list = [
        (
        protein,
        ligand_pdbqt_path,
        config_file,
        os.path.join(
            "/Users/yalilyu/Desktop/code/DrugSynthMC2/docking/ar/ar_structures",
            f"{protein}_clusters",
            f"{protein}_clusters_{i}.pdbqt"
        ),
        i
        )
        for i in range(1, num_files +1)
    ]

    with Pool() as pool:
        results = pool.map(run_vina, args_list)
    
    with open(log_path, 'a') as log:
        for result in results:
            log.write(result)
    
    end_time = time.time()
    print("Docking time used", end_time - start_time)
    
    affinities = []
    in_table = False

    with open(log_path, 'r') as f:
        for line in f:
            line_strip = line.strip()

            # Detect start of a table
            if line_strip.startswith('--- Affinity Table'):
                in_table = True
                continue

            # Detect separator line that precedes data
            if in_table and line_strip.startswith('-----+'):
                continue

            # Exit table when encountering a new table
            if in_table and line_strip.startswith('--- Affinity Table'):
                in_table = False
                continue

            # Collect affinity values
            if in_table:
                parts = line_strip.split()
                if len(parts) >= 2:
                    try:
                        affinity_value = float(parts[1])
                        affinities.append(affinity_value)
                    except ValueError:
                        pass

    return min(affinities) if affinities else None


# def vina_docking(protein, ligand_pdbqt_path, config_file, num_files):
#     """
#     Convert ligand smiles string to pdbqt, then batch Docking the ligand to cluster structures of a specified protein, write affinties table to log, then extract the most favourable affinity value.
#     """
#     start_time = time.time()
    
#     with open(f'/Users/yalilyu/Desktop/python_test_docking/{protein}.txt', 'w') as file:
#     # erase log contents
#         pass

#     for i in range(1, num_files + 1):
#         protein_file = os.path.join(
#             "/Users/yalilyu/Desktop/python_test_docking/",
#             f"{protein}_clusters",
#             f"{protein}_clusters_{i}.pdbqt"
#         )

#         # Run Vina and capture output
#         command = [
#         "vina",
#         "--receptor", protein_file,
#         "--ligand", ligand_pdbqt_path,
#         "--config", config_file
#         ]

#         result = subprocess.run(command, capture_output=True, text=True, check=True)
#         stdout = result.stdout

#         # Extract only the affinity table
#         affinities_table = extract_affinity_table_from_stdout(stdout)
    
#         log_path = f'/Users/yalilyu/Desktop/python_test_docking/{protein}.txt'
    
#         if not os.path.exists(log_path):
#             with open(log_path, 'w') as f:
#                 pass  # create an empty file
#             print(f"Created empty file: {log_path}")

#     # Append only the table into log
#         with open(log_path, 'a') as log:
#             log.write(f"\n--- Affinity Table for {protein} Cluster {i} ---\n")
#             log.write(affinities_table)
#             log.write('\n\n')

#         print(f"{protein} cluster {i} results appended to log.")
    
#     end_time = time.time()
#     print("Docking time used", end_time - start_time)

#     affinities = []
#     in_table = False

#     with open(log_path, 'r') as f:
#         for line in f:
#             line_strip = line.strip()

#             # Detect start of a table
#             if line_strip.startswith('--- Affinity Table'):
#                 in_table = True
#                 continue

#             # Detect separator line that precedes data
#             if in_table and line_strip.startswith('-----+'):
#                 continue

#             # Exit table when encountering a new table
#             if in_table and line_strip.startswith('--- Affinity Table'):
#                 in_table = False
#                 continue

#             # Collect affinity values
#             if in_table:
#                 parts = line_strip.split()
#                 if len(parts) >= 2:
#                     try:
#                         affinity_value = float(parts[1])
#                         affinities.append(affinity_value)
#                     except ValueError:
#                         pass

#     return min(affinities) if affinities else None

def docking_score(protein, smile_string, config_file, num_files):
    """
    Docking the generated ligand onto the protein and return the lipinski score as a component of the reward, and erase log
    """
    ligand_pdbqt_path, ligand_pdbqt_content = prepare_ligand(smile_string)
    best_affinity = vina_multiprocessing(protein, ligand_pdbqt_path, config_file, num_files)
    affinity_score = 1 - 1 / (1 + math.exp(-(best_affinity +5)))
    
    return best_affinity, affinity_score

# def docking_score(protein, smile_string, config_file, num_files):
#     """
#     Docking the generated ligand onto the protein and return the lipinski score as a component of the reward, and erase log
#     """
#     ligand_pdbqt_path, ligand_pdbqt_content = prepare_ligand(smile_string)
#     best_affinity = vina_docking(protein, ligand_pdbqt_path, config_file, num_files)
#     affinity_score = 1 - 1 / (1 + math.exp(-(best_affinity +5)))
    
#     return best_affinity, affinity_score

