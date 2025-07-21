from rdkit import Chem
from rdkit.Chem import BRICS, Descriptors, AllChem
from rdkit.Chem.rdmolfiles import MolFromSmiles, MolToSmiles
from rdkit.Geometry import Point3D
import re
import numpy as np
from rdkit.Chem import rdChemReactions as Reactions



# def BondIdxFromAtomIdx(mol, atom1_idx, atom2_idx):
#     for bond in mol.GetBonds():
#         a1 = bond.GetBeginAtomIdx()
#         a2 = bond.GetEndAtomIdx()
#         if (a1 == atom1_idx and a2 == atom2_idx) or \
#         (a2 == atom1_idx and a1 == atom2_idx):
#             BRICS_bond = int((bond.GetIdx()))
#             return BRICS_bond
#         else:
#             continue
#     return None


# def GetBRICSBonds(mol):
#     BRICS_bonds = []
#     atomPairs_and_envs = list(BRICS.FindBRICSBonds(mol))

#     for atomPair_and_env in atomPairs_and_envs:
#         atomPair = atomPair_and_env[0]
#         atom1_idx = atomPair[0]
#         atom2_idx = atomPair[1]
#         BRICS_bond = BondIdxFromAtomIdx(mol,atom1_idx,atom2_idx)
#         BRICS_bonds.append(BRICS_bond)
        
#     return BRICS_bonds
    
    

def getSingleFrag(frags):
    """
    After molecule decomposition, retain the frag with a larger molecular weight
    """
    frags = list(frags)
    retained_frag_wt = Descriptors.ExactMolWt(frags[0])
    retained_frag = frags[0]
    for frag in frags:
        wt = Chem.Descriptors.ExactMolWt(frag)
        if wt > retained_frag_wt:
            retained_frag = frag
    
    return retained_frag



def getRetainedFrags(smiles):
    """
    Get all retained (larger) fragments, along with their smiles, from different ways (defined by BRICS) of molecule decomposition
    """
    # get number of bonds that can be broken
    mol = MolFromSmiles(smiles)
    bond_list = list(BRICS.FindBRICSBonds(mol))
    num_BRICSBonds = len(bond_list)

    retained_frags = []
    for  i in range(num_BRICSBonds):
        retained_frag = []
        # by specifying the bond to be broken in BRICS.BreakBRICSBonds, we can break one bond at a time
        bond = bond_list[i]
        broken_mol = BRICS.BreakBRICSBonds(mol,bonds=[bond]) 
        atom_idx1, atom_idx2 = bond[0]
        evn1, evn2 = bond[1]
        frag_mol = Chem.GetMolFrags(broken_mol,asMols=True)
        retained_mol = getSingleFrag(frag_mol)
        retained_frag_smile = MolToSmiles(retained_mol, kekuleSmiles = True)

        frag_env = re.match(r"^\[(\d+)\*",retained_frag_smile).group(1)
        retained_mol._env = frag_env
        if frag_env == evn1:
            # retained_mol._atomidx = atom_idx1
            retained_mol._compenv = evn2
            # retained_mol._compatomidx = atom_idx2
        else:
            # retained_mol._atomidx = atom_idx2
            retained_mol._compenv = evn1
            # retained_mol._compatomidx = atom_idx1

        retained_frag.append(retained_mol)
        retained_frag.append(retained_frag_smile)
        retained_frags.append(retained_frag)

    return retained_frags

