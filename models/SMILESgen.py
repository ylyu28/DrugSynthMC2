from dataclasses import dataclass
from typing import List, ClassVar, Dict, Set, Tuple, Any
import json
import copy 
from rdkit import Chem
from rdkit.Chem import Lipinski
from tools.NNreader import Prediction, Model
import math
from docking.docking import docking_score
import re

# S(=O)(=O) -> U
# S(=O) -> M
# Cl -> L
# C(F)(F)(F) -> W
ATOMS = ['C', 'O', 'N', 'F', 'S', 'P', 'U', 'M', 'L', 'W']

NUMBERS = ['1', '2', '3', '4', '5', '6', '7', '8', '9']

NGRAM_LEN = 3 # at least 1

PRUNING_MOVE_TR = 0.001 # frequency under which we prune the move, deemed as too rare

MAX_RING_LEN = 7

MIN_RING_LEN = 5

HEURISTIC_MODE = 'ngram' # 'ngram' or 'neural'

# @dataclass
# class NA:
#     nextAtom: Dict[str, float]
    
#     def __init__(self, data: Dict[str, float]):
#         self.nextAtom = data
    
#     def __repr__(self):
#         return f"NA(nextAtom={self.nextAtom})"

legal_bonds: frozenset[Tuple[str,str,int]] = frozenset({
    ('C', 'C', 1),
    ('C', 'O', 1),
    ('C', 'N', 1),
    ('C', 'F', 1),
    ('C', 'S', 1),
    ('C', 'P', 1),
    ('C', 'B', 1),
    ('C', 'c', 1),
    ('C', 'I', 1),
    ('O', 'N', 1),
    ('N', 'N', 1),
    ('N', 'S', 1),
    ('N', 'P', 1),
    ('P', 'S', 1),
    ('C', 'C', 2),
    ('C', 'O', 2),
    ('C', 'N', 2),
    ('C', 'S', 2),
    ('O', 'N', 2),
    ('O', 'S', 2),
    ('N', 'N', 2),
    ('N', 'P', 2),
    ('S', 'P', 2),
})


prior = Model.load_model("Neural/SMILESexplicit_shortcuts")

# ngrams: Dict[str, NA] ={}
ngrams_path = "ngrams/fda_ngrams_shortcuts_cycles.json"

with open(ngrams_path, 'r') as f:
    ngrams = json.load(f)


@dataclass(eq=True)
class Move:
    atom: str
    doubleLink: bool
    nesting: bool
    closeNesting: bool
    cycle: int

    def __hash__(self):
        return hash((self.atom, self.doubleLink, self.nesting, self.closeNesting, self.cycle))

class State:
    CONSIDER_NON_TERM: bool = False
    BEST_POSSIBLE_SCORE: float = 2.924 # at which lipinskiness fulfilled and kd < -9.5

    def __init__(self):
        self.SMILE = ['C']
        self.nestingOpenAtom = ['C']
        self.nestingOpenCovalence = [4]
        self.nestingCycleToClose = [0]
        self.openCycles = []
        self.openCyCosts = []
        self.closedCycles = 0
        self.seq = []
        self.open_nesting_ASAP = False
        self.finish_ASAP = False
        self.reached_best_score = False
        self.target_OtoC_ratio = -1.0
        self.target_NtoC_ratio = -1.0
        self.storedPrior = Prediction(label=[],confidence=[])
    
    @classmethod
    def new(cls):
        return cls()
    
    def play(self, m: Move):
        addition = ' '
        if m.doubleLink and m.cycle == 0: # when m.cycle != 0, we are closing a ring which started with double bond, (i.e. with TWO sequential moves 'number' + '=', since ring opening cannot be accompanied by doulbeLink), therefore no need to add '=' to the smile. 
            addition = '='

            if self.SMILE[len(self.SMILE) - 1] == '(':
                self.nestingOpenCovalence[-2] -= 1
            else:
                self.nestingOpenCovalence[-1] -= 1

        if m.nesting:
            addition = '('
            self.nestingOpenCovalence[-1] -= 1
            self.nestingOpenCovalence.append(1)
            self.nestingOpenAtom.append('X')

            if self.open_nesting_ASAP:
                self.open_nesting_ASAP = False
                self.nestingCycleToClose.append(self.openCycles[-1])
            else:
                self.nestingCycleToClose.append(0)
        
        if m.closeNesting:
            addition = ')'
            self.nestingOpenCovalence.pop()
            self.nestingOpenAtom.pop()
            self.nestingCycleToClose.pop()
        
        if m.cycle != 0:
            addition = str(m.cycle)
            
            if m.cycle in self.openCycles:
                self.closedCycles += 1
                self.openCycles.remove(m.cycle)

                for i in range(len(self.nestingCycleToClose)):
                    if self.nestingCycleToClose[i] == m.cycle:
                        self.nestingCycleToClose[i] = 0
            else:
                self.open_nesting_ASAP = True
                self.openCycles.append(m.cycle)

            self.nestingOpenCovalence[-1] -= 1
            
            if m.doubleLink:
                self.nestingOpenCovalence[-1] -= 1
        
        if m.atom != ' ':
            addition = m.atom

            new_covalence = 0
            if m.atom == 'C':
                new_covalence = 3
            if m.atom == 'N':
                new_covalence = 2
            if m.atom == 'O':
                new_covalence = 1
            if m.atom == 'F':
                new_covalence = 0
            if m.atom == 'S':
                new_covalence = 1
            if m.atom == 'U':
                new_covalence = 1
            if m.atom == 'M':
                new_covalence = 1
            if m.atom == 'L':
                new_covalence = 0
            if m.atom == 'W': 
                new_covalence = 0
            
            for i in range(len(self.SMILE)):
                last_char = self.SMILE[-1-i]

                if last_char == '=':
                    new_covalence -= 1
                else:
                    if last_char != '(':
                        if last_char != ')':
                            break
                    
            self.nestingOpenCovalence[-1] = new_covalence
            self.nestingOpenAtom[-1] = m.atom
        
        self.SMILE.append(addition)

        if len(self.SMILE) > 60: # to be discussed
            self.finish_ASAP = True
        
        self.seq.append(m)

                                                
    def legal_moves(self) -> list:
        """
        Defining 6 legal moves: 
        1. closing a ring
        2. closing a nesting
        3. opening a nesting
        4. adding an atom
        5. adding a double bond
        6. opening a ring
        """
        legal_moves = []
        last_char = self.SMILE[-1]
        second_last = 'x'
        if len(self.SMILE) > 1:
            second_last = self.SMILE[-2]
        
        # checking if there is an open nesting level that could continue afterwards. 
        # If there isn't, end with a fluorine, or a double bond and a sulfur or oxygen, 
        # while having an open cycle is prohibited.
        could_prevent_cycle_completion = True
        for i in range(len(self.nestingOpenCovalence)-1):
            # element in nestingCycleToClose = 0 when the cycle is closed (ie when the cycle number shows its second presence in the smile string), this 0 then popped when nesting closed
            if self.nestingCycleToClose[len(self.nestingOpenCovalence) -1-i] != 0:  # While there is at least one open cycle, there is a possibility that the cycle might be prevented from closing
                break
            if self.nestingOpenCovalence[-2-i] != 0: # when the cycle is already closed, if the last open covalence atom has at least one available bond, the molecule has capacity to continue growing
                could_prevent_cycle_completion = False
        
        can_play_end_cycle = False
        must_close_cycle = False

        # Legal move 1: ring closure (only possible when at least one ring is open)
        if len(self.openCycles) >= 1 : 
        # if len(self.openCycles) >= 1: # if at least one ring is open
            cycle_len, _, proper_atoms, bond_cost = self.backtrackCycle(self.SMILE.copy(), str(self.openCycles[-1]))
            # bond_cost = valence needed to close the ring
            # if more than one ring is open, we need to make sure the atom at which the current ring is growing has capacity to close current ring + grow previous ring
            if not self.finish_ASAP or (len(self.openCycles) == 1 or self.nestingOpenCovalence[-1] >= 1 + bond_cost): # nestingOpenCovalence[-1] is the atom at which the cycle is growing
            
                if last_char != '=' and not (last_char == '(' and second_last == '='): # only checking last_char and second_last if last_char != '='
                    # last_char != '=': 
                    # "AIzynthfinder forbids "=(", only "(=" is allowed."
                    
                    # "Closing a cycle is prohibited if it blocks everything afterwards."
                    if not (
                        len(self.openCycles) > 1 and could_prevent_cycle_completion and self.nestingOpenCovalence[-1] <= bond_cost
                        ) and (self.nestingOpenCovalence[-1] >= 1 + bond_cost - 1):

                        if cycle_len >= MIN_RING_LEN and proper_atoms >= 4:
                            mv = Move(atom = ' ', doubleLink = False, nesting = False, closeNesting = False, cycle = self.openCycles[-1])
                            if bond_cost == 2:
                                mv.doubleLink = True
                            legal_moves.append(mv)
                            can_play_end_cycle = True

                        if cycle_len >= MAX_RING_LEN:
                            must_close_cycle = True

            # if below is true, the only legal_move considered is ring closure
            if must_close_cycle and can_play_end_cycle and (not could_prevent_cycle_completion or self.nestingOpenCovalence[-1] >= 1 + bond_cost) and self.nestingOpenCovalence[-1] >= bond_cost:
                return legal_moves
            
        # Legal move 2: close a nesting, prohibited when right after a nesting or a double link, or while a cycle is not closed yet
        if len(self.nestingOpenCovalence) != 1 and self.nestingCycleToClose[-1] == 0 and not self.open_nesting_ASAP:
            if last_char != '=' and last_char != '(' and not (could_prevent_cycle_completion and len(self.openCycles) >0):
                mv = Move(atom = ' ', doubleLink = False, nesting = False, closeNesting = True, cycle = 0)
                legal_moves.append(mv)
        
        # again: immediate RING closure as nesting closure is not considered as a legal move when self.nestingCycleToClose[-1] != 0
        if must_close_cycle and self.nestingCycleToClose[-1] != 0:
            return legal_moves
        
        # Legal move 3: open a nesting which is the only legal move when the 'open_nesting_ASAP' flag == True
        # "AIzynthfinder forbids "=(", only "(=" is allowed."
        if self.nestingOpenCovalence[-1] >= 1 or last_char == '(':
                # opening of a nesting, prohibited in finish ASAP
            if last_char != '(' and last_char != '=' and (not self.finish_ASAP or (must_close_cycle and not can_play_end_cycle)):
                mv = Move(atom = ' ', doubleLink = False, nesting = True, closeNesting = False, cycle = 0)
                legal_moves.append(mv)
                
            if self.open_nesting_ASAP and (last_char in ATOMS):
                # if open_nesting_ASAP, legal_move only includes open nesting
                mv = Move(atom = ' ', doubleLink = False, nesting = True, closeNesting = False, cycle = 0)
                return [mv]
            
        # Legal move 4: adding an atom
        if self.nestingOpenCovalence[-1] >= 1 or last_char == '(':
            if not self.finish_ASAP or last_char in ('=', '(') or len(self.openCycles) > 0 or (self.open_nesting_ASAP and (last_char not in ATOMS)): # "(self.openCycles.len() > 0 && !can_play_end_cycle) removed the clause on the end of cycle to avoid having too many forced cycles of size 5."
                
                for i in ['C', 'O', 'F', 'N', 'S', 'U', 'M', 'L', 'W']:
                    # "We also sauté if it's F and we're at the top level and there are still loops to close."
                    # "We prohibit S and O if we are in finish ASAP unless it completes an = and does not prevent finishing the cycles, or if it could block the end of the cycles."
                    
                    # conditions for F, L -> Cl, W -> C(F)(F)(F) that would prohibit addition of the three atoms # single valence
                    condition1 = (
                        self.nestingCycleToClose[-1] != 0 or # a ring is opened and not closed yet
                        (len(self.openCycles) > 0 and could_prevent_cycle_completion) or
                        self.open_nesting_ASAP or
                        last_char == '=' or
                        (last_char == '(' and second_last == '=')
                    )
                    # conditions for O, U -> S(=O)(=O), M -> S(=O) # double valence
                    condition2 = (
                        (self.finish_ASAP and last_char != '=' and not could_prevent_cycle_completion) or 
                        (len(self.openCycles) > 0 and could_prevent_cycle_completion) or
                        (last_char == '=' and self.nestingCycleToClose[-1] != 0) or # ie the atom to be added needs capacity to grow the cycle
                        (last_char == '=' and self.open_nesting_ASAP) # the atom to be added needs capacity to open nesting
                    )

                    if not ((i in ('F','L','W') and condition1) or (i in ('O','U','M') and condition2)):
                        
                        prev_atom = self.nestingOpenAtom[-1] # previous atom that the new atom is to bond to 
                        n = 0
                        while prev_atom == 'X': # placeholder atom when nesting opened
                            n += 1
                            prev_atom = self.nestingOpenAtom[len(self.nestingOpenCovalence) -1-n]
                            
                        bondType = 1
                        if last_char == '=':
                            bondType = 2

                        if ((i, prev_atom, bondType) in legal_bonds) or ((prev_atom, i, bondType) in legal_bonds) or i == 'U' or i == 'M' or i == 'L' or i == 'W': 
                            mv = Move(atom = i, doubleLink = False, nesting = False, closeNesting = False, cycle = 0)
                            legal_moves.append(mv)

        # Legal move 5: adding a double bond
        # Double link and cycle opening are prohibited in finish ASAP."
        if (self.nestingOpenCovalence[-1] >= 2 or (last_char == '(' and self.nestingOpenCovalence[-2] >= 1)) and not self.finish_ASAP:
             # Double link is prohibited after a cycle starts.
            if last_char != '=':
                mv = Move(atom = ' ', doubleLink = True, nesting = False, closeNesting = False, cycle = 0)
                legal_moves.append(mv)

        # Legal move 6: opening a cycle
            # "Opening of a cycle, maximum 9."
            # if len(self.openCycles) + self.closedCycles < 9 and last_char != '(' and (last_char not in NUMBERS): # Initially, it stopped at '='
            if len(self.openCycles) + self.closedCycles < 9 and last_char != '(' and (last_char not in NUMBERS) and len(self.openCycles) <=2: # May 29, constraining the software from opening a third ring two already open to avoid weirdly concatenated rings
                mv = Move(atom = ' ', doubleLink = False, nesting = False, closeNesting = False, cycle = len(self.openCycles) + self.closedCycles + 1)
                legal_moves.append(mv)

        if PRUNING_MOVE_TR >= 0.0:
            legal_moves_copy = legal_moves.copy()
            vecbackup = legal_moves.copy()
            legal_moves.clear()

            for mv in legal_moves_copy:
                heuri = self.heuristic(mv)
                if heuri != 0.0 and math.exp(heuri) > PRUNING_MOVE_TR:
                    legal_moves.append(mv)

            if len(legal_moves) == 0 and len(vecbackup) > 0:
                legal_moves = vecbackup
            
        return legal_moves
                      


    def terminal(self) -> bool:
        return len(self.legal_moves()) == 0
    

    
    def score(self) -> float:
        try:
            kd, aff_sc = docking_score('AR', self.smile_to_smile(self.SMILE), "./docking/ar/ar_box.txt",1)
        except:
            return 1000, -2000
        else:
            sc = aff_sc + self.lipinskiness()
            if sc >= self.BEST_POSSIBLE_SCORE:
                self.reached_best_score = True
            return kd, sc
        
             
    def backtrackCycle(self, SMILE: list, last_open_cycle: str) -> tuple:
        cycle_length = 0
        left_current_nesting_level = 0
        right_current_nesting_level = 0
        left_chain = []
        right_chain = []
        left_chain_id = []
        right_chain_id = []
        # would_be_made_rigid = []

        cycle_encountered = False
        improper_atoms = 0
        active_subcycles = []

        bond_cost = 1
        s = ""
        for c in SMILE:
            s += c
        
        for i in range(len(SMILE)):
            indice = len(SMILE) - 1 - i

            if SMILE[indice] == last_open_cycle:
                cycle_encountered = True # ring opened

                if SMILE[indice-1] == '=':
                    bond_cost = 2 # the current nesting atom is to be connected to a double bond to close the ring
            
            if SMILE[indice] in NUMBERS:
                if SMILE[indice] in active_subcycles:
                    active_subcycles.pop()
                else:
                    active_subcycles.append(SMILE[indice])

            if SMILE[indice] == ')':
                right_current_nesting_level -= 1
                if cycle_encountered:
                    left_current_nesting_level -= 1
            
            if SMILE[indice] == '(':
                right_current_nesting_level += 1
                if right_current_nesting_level > 0:
                    right_current_nesting_level = 0
                if cycle_encountered:
                    left_current_nesting_level += 1
                    if left_current_nesting_level > 0: 
                        left_current_nesting_level = 0
            
            if SMILE[indice] in ATOMS:
                if left_current_nesting_level == 0 and right_current_nesting_level == 0 and cycle_encountered:
                    chain = right_chain.copy()
                    chain_id = right_chain_id.copy()
                    for i in range(len(left_chain)):
                        chain.append(left_chain[-1-i])
                        chain_id.append(left_chain_id[-1-i])

                    aromatic = True
                    last_equal = 1
                    if len(chain) != 0 and chain[0] == '=':
                        last_equal = 0
                    
                    double_count = 0
                    atom_count = 1
                    for i in range(len(chain)):
                        c = chain[i]
                        if c == '=':
                            # would_be_made_rigid.pop()

                            double_count += 1
                            if last_equal < 2:
                                aromatic = False
                            
                            if last_equal > 2:
                                aromatic = False

                            last_equal = 0
                        if c in ATOMS:
                            last_equal += 1
                            atom_count += 1
                            # would_be_made_rigid.append(chain_id[i])
                    
                    if double_count * 2 < atom_count -1:
                        aromatic = False
                    
                    # if not aromatic:
                        # would_be_made_rigid = []
                    
                    # return (cycle_length+1, aromatic, would_be_made_rigid, cycle_length - improper_atoms, bond_cost)
                    return (cycle_length+1, aromatic, cycle_length - improper_atoms, bond_cost)
                
                if left_current_nesting_level == 0 and cycle_encountered:
                    left_chain.append(SMILE[indice])
                    left_chain_id.append(indice)
                    cycle_length += 1
                    if len(active_subcycles) > 0:
                        improper_atoms += 1
                
                if right_current_nesting_level == 0:
                    right_chain.append(SMILE[indice])
                    right_chain_id.append(indice)
                    cycle_length += 1
                    if len(active_subcycles) > 0:
                        improper_atoms += 1
            
            if SMILE[indice] == '=':

                if left_current_nesting_level == 0 and cycle_encountered:
                    left_chain.append(SMILE[indice])
                    left_chain_id.append(indice)
                
                if right_current_nesting_level == 0:
                    right_chain.append(SMILE[indice])
                    right_chain_id.append(indice)
            
        print("Error, backtrack cycle should not return here")
        print(s)
        # would_be_made_rigid = []
        return (cycle_length, False, cycle_length - improper_atoms, bond_cost)
        
    def make_from_string(s:str):
        """
            
        """
        st = State.new()
        for c in s:
            m = Move(atom = ' ', doubleLink = False, nesting = False, closeNesting = False, cycle = 0)

            if c == '(':
                m.nesting = True
            if c == ')':
                m.closeNesting = True
            if c == '=':
                m.doubleLink = True
            if c in NUMBERS:
                m.cycle = int(c)
            if c in ATOMS:
                m.atom = c
                
            st.play(m)
            
        return st
    
    def smile_to_smile(self, SMILE: list) -> str:
        """
        Convert SMILE list to real smile string, whilst converting shortcut representations.
        """
       
        smile = ""
        for i in SMILE:
            if i == 'U':
                smile +='S(=O)(=O)'
            elif i == 'M':
                smile += 'S(=O)'
            elif i == 'L':
                smile += 'Cl'
            elif i == 'W':
                smile += 'C(F)(F)(F)'
            else:
                smile += i
        
        smile = re.sub(r'\((\d)\)',r'\1',smile) # fixing the (ring_number) pattern
                
        return smile

    # def get_mol_counts(self, SMILE: list) -> tuple:
    #     """
    #     Input self.SMILE and get total number of atoms in the molecule
    #     """

    #     smile = self.smile_to_smile(SMILE)
    #     mol = Chem.MolFromSmiles(smile)
    #     mol = Chem.rdmolops.AddHs(mol)
    #     n_atoms = Chem.rdMolDescriptors.CalcNumAtoms(mol)
    #     return n_atoms
    
        # num_h = 0
        # for a in mol.GetAtoms():
        #     num_h += Chem.rdchem.GetTotalNumHs(a)
        # num_heavy = Chem.rdMolDescriptors.CalcNumHeavyAtoms(mol)
        # return (num_heavy, num_h)

    # def ratioH(self, SMILE: list) -> float:
    #     """
    #     Input self.SMILE and get the ratio of number (heavy atoms) / number (hydrogens)

    #     """

    #     (num_heavy, num_h) = self.get_mol_counts(SMILE)
    #     ratioH = num_heavy / num_h
    #     return ratioH


    def lipinskiness(self) -> float:
        """
        Calculate Lipinskiness score for a generated molecule, which is the total score of 9 sub-scores,
        best possible Lipinskiness score is 2
        """
        smile = self.smile_to_smile(self.SMILE)
        mol = Chem.MolFromSmiles(smile)

        if mol is None:
            return -2.0
        
        mol = Chem.rdmolops.AddHs(mol)

        n_hbond_donor = Lipinski.NumHDonors(mol)
        n_hbond_acceptor = Lipinski.NumHAcceptors(mol)
        n_atoms = Chem.rdMolDescriptors.CalcNumAtoms(mol)
        molecular_weight = Chem.rdMolDescriptors.CalcExactMolWt(mol)
        nitro_count = 0
        carbon_count = 0
        oxygen_count = 0
        n_rigid_bonds = 0
        n_rings = 0
        for i in self.SMILE:
            if i == 'N':
                nitro_count += 1
            if i == 'C':
                carbon_count += 1
            if i == 'O' or i == 'U' or i == 'M':
                oxygen_count += 1
            if i == '=':
                n_rigid_bonds += 1
            if i in NUMBERS:
                n_rings += 0.5

        # Calculate lipinski score
        lipinski_sc = 0
        # 1st rule: < 5 hydrogen bond donors
        lipinski_sc += -max(n_hbond_donor - 5.0, 0.0)/5.0 # highest score (when n_hbond_donor < 5) is 0.0
        # 2nd rule: < 10 hydrogen bond acceptors
        lipinski_sc += -max(n_hbond_acceptor - 10.0, 0.0)/10.0 # highest score (when n_hbond_accptor < 10) is 0.0
        # 3rd rule: < 800 molecular weight
        lipinski_sc += -max(molecular_weight-800.0, 0.0)/800.0 # highest score (when molecular_weight < 800) is 0.0
        # 4th rule: > 30 number of atoms
        lipinski_sc += min(n_atoms, 30.0)/30.0 # highest score = 1 (when n_atoms > 30)
        # 5th rule: < 70 number of atoms
        lipinski_sc += -max(n_atoms - 70.0, 0.0)/70.0 # highest score = 0 (when n_atoms < 70)
        # 6th rule: max score = 1 (at least one rigid bond for every 6 atoms)
        lipinski_sc += min(n_rigid_bonds, n_atoms/6.0)/(n_atoms/6.0)
        # 7th rule: max score = 0 (when number of rings < 5)
        lipinski_sc += -max(n_rings - 5.0, 0.0)/5.0
        # 8th rule: self.target_NtoC_ratio - NtoC_ratio).abs() < 0.1
        if self.target_NtoC_ratio >= 0.0:
            NtoC_ratio = nitro_count/carbon_count
            lipinski_sc += -max((self.target_NtoC_ratio - NtoC_ratio).abs() - 0.1, 0.0) # highest score = 0 (when self.target_NtoC_ratio - NtoC_ratio).abs() < 0.1)
        # 9th rule: self.target_OtoC_ratio - OtoC_ratio).abs() < 0.1
        if self.target_OtoC_ratio >= 0.0:
            OtoC_ratio = oxygen_count/carbon_count
            lipinski_sc += -max((self.target_OtoC_ratio - OtoC_ratio).abs() - 0.1, 0.0) # highest score = 0 (when self.target_OtoC_ratio - OtoC_ratio).abs() < 0.1)
        # yali: new rule re ring not closed yet issue
        # if len(self.openCycles) != 0:
        #     lipinski_sc += -2
        
        return lipinski_sc 
    
    
    
    def heuristic(self, m: Move) -> float:
        mv = 'a'
        if m.doubleLink:
            mv = '='
        if m.nesting:
            mv = '('
        if m.closeNesting:
            mv = ')'
        if m.cycle !=0:
            mv = 'X'
        if m.atom != ' ':
            mv = m.atom

        if HEURISTIC_MODE == 'neural':
            SMILEstring = []
            for c in self.SMILE:
                SMILEstring.append(c)
            
            SMILEstr = []
            for c in SMILEstring:
                SMILEstr.append(c)

            if len(self.storedPrior.label) == 0:
                self.storedPrior = prior.predict(SMILEstr)
        
            pred = self.stroedPrior.copy()

            val = ["\n", "&", "C", "(", ")", "1", "=", "2", "O", "N", "3", "F", "[C@@H]", "#", "S", "L", "[O-]", "[C@H]", "[NH+]", "[C@]", "Br", "/", "[NH3+]", "W", "4", "[NH2+]", "U", "[C@@]", "[N+]", "\\", "M", "[S@]", "5", "[N-]", "[S@@]", "[S-]", "6", "7", "I", "P", "[OH+]", "[NH-]", "[P@@H]", "[P@@]", "[PH2]", "[P@]", "[P+]", "[S+]", "[O+]", "[CH2-]", "[CH-]", "[SH+]", "[PH+]", "[PH]", "8", "[S@@+]"]

            if mv == 'X':
                sum_conf = 0
                for i in range(len(pred.label)):
                    try:
                        float(val[pred.label[i]])
                        sum_conf += float(pred.confidence[i])
                    except ValueError:
                        pass
                return math.log(sum_conf)
            
            if mv == m.atom:
                sum_conf = 0
                for i in range(len(pred.label)):
                    if val[pred.label[i]] == str(mv) or val[pred.label[i]] == mv.lower():
                        sum_conf += float(pred.confidence[i])
                return math.log(sum_conf)
            
            for i in range(len(pred.label)):
                if val[pred.label[i]] == str(mv):
                    return math.log(pred.condifence[i])
                
        if HEURISTIC_MODE == "ngram":
            s = ""
            if mv == 'X' and (m.cycle in self.openCycles):
                (si, _, _, _) = self.backtrackCycle(self.SMILE.copy(),str(m.cycle))
                if si == 5:
                    return math.log(0.595046/(0.595046+1.958514+0.053870))
                if si == 6:
                    return math.log(1.958514/(1.958514+0.053870))
                if si == 7:
                    return math.log(1.0)
            
            for i in range(NGRAM_LEN):
                if len(self.SMILE) > NGRAM_LEN -i -1:
                    if self.SMILE[-(NGRAM_LEN-i)] in NUMBERS:
                        s += 'X'
                    else:
                        s += self.SMILE[-(NGRAM_LEN-i)]

            key = ngrams.get(s)
            if key is None:
                return 0.0
            
            ret = key.get(mv)
            if ret is None:
                return 0.0
            
            return math.log(ret)
        return 0.0
    

    def soothedScore(self) -> float:
        return self.score()


            
            







                    




            


        


        
