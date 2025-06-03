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
from tools.resultSaver import writeline

def prepare_ligand(smile_string):
    print("start preparing ligand")
    """
    molscrub used to assign protonation states and generate ligand conformers
    writer functions in meeko to save the prepared ligand molecules into PDBQT files
    """
    pdbqt_dir = "./result/pdbqt"
    os.makedirs(pdbqt_dir, exist_ok=True)
    
    pdbqt_path = f"./result/pdbqt/{smile_string}.pdbqt"
    ligand_sdf_path = f"./result/pdbqt/ligand.sdf"
    
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
    print("sdf ready")
    subprocess.run([
        "mk_prepare_ligand.py",
        "-i", ligand_sdf_path,
        "-o", pdbqt_path
    ], check=True)

    with open(pdbqt_path, 'r') as f:
        ligand_content = f.read()
    print("ligand prepared")
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
    "--config", config_file,
    "--num_modes","1"
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
    log_path = f'./docking/{protein}.txt'

    with open(log_path, 'w') as f:
        pass  # create an empty file
    
    args_list = [
        (
        protein,
        ligand_pdbqt_path,
        config_file,
        os.path.join(
            "./docking/ar/ar_structures",
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
    time_taken = time.time() -start_time
    print(f"min affinity is {min(affinities)}, multiprocessing time taken {time_taken}")
    return min(affinities) if affinities else None

# def docking_kd(protein, smile_string, config_file, num_files):
#     ligand_pdbqt_path, ligand_pdbqt_content = prepare_ligand(smile_string)
#     best_affinity = vina_multiprocessing(protein, ligand_pdbqt_path, config_file, num_files)
#     return best_affinity

def docking_score(protein, smile_string, config_file, num_files):
    """
    Docking the generated ligand onto the protein and return the lipinski score as a component of the reward, and erase log
    """
    # start_time = time.time()
    ligand_pdbqt_path, ligand_pdbqt_content = prepare_ligand(smile_string)
    best_affinity = vina_multiprocessing(protein, ligand_pdbqt_path, config_file, num_files)
    affinity_score = 1 - 1 / (1 + math.exp(-(best_affinity +7)))
    # writeline(str(time.time() - start_time) + " "+ smile_string + " " + str(best_affinity)+"\n","affinityMonitor")
    return best_affinity, affinity_score

