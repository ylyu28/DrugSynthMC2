from models import SMILESgen_cyanide
from methods import fragNMCS
import time
import fragments.Frags as Frags
import rdkit.Chem as Chem

if __name__ == '__main__':
    molGenState = SMILESgen_cyanide.State.new()

    print("Starting NMCS search...")
    frags = Frags.getRetainedFrags('CC(CN1C=CC(=N1)C2=CC(=C(C=C2)C#N)Cl)NC(=O)C3=NNC(=C3)C(C)O') 
    frag_mols = [frag[0] for frag in frags]
    frag_smiles = [Chem.MolToSmiles(mol, kekuleSmiles=True) for mol in frag_mols]
    print(frag_smiles)
    
    for frag_mol in frag_mols:
        num_molecules = 0
        while num_molecules < 100:
            st = fragNMCS.launch_nmcs('ar',frag_mol, molGenState, level=2, heuristic_w=1.0, verbose=True, timeout=0.0, register_name='test')
            num_molecules +=1


