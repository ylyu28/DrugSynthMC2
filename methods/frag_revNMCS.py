from models.SMILESgen_cyanide import State
from tools.calc import softmaxChoice
from tools.resultSaver import writeline
import time
import random
import copy
import math


class frag_revNMCS:
    def __init__(self):
        self.best_yet = -1000.0
        self.timeout = -1.0
        self.registerName= ""
        self.start_time = time.perf_counter()
    
    @classmethod
    def new(cls):
        return cls()
    

    def playout(self, st:State, heuristic_w):
        while not st.terminal():
            moves = st.legal_moves()
            if len(moves) == 0:
                break
            
            i = int(len(moves) * random.random())
            if heuristic_w != 0.0:
                weights = [heuristic_w * st.heuristic(m) for m in moves]
                i = softmaxChoice(weights)
            
            mv = moves[i]
            st.play(mv)
        
        return st
    

    
    def nmcs(self, protein, frag_mol, st: State, level, heuristic_w, verbose: bool, initial_level):
        best_state = copy.deepcopy(st)
        best_state_lip = -1000
        best_state_kd = 1000

        while not st.terminal(): # runs until the state is terminal (no legal moves are left)
            moves = st.legal_moves()
            if len(moves) == 0:
                break
            for mv in moves:
                if (time.time() - self.start_time) > self.timeout and self.timeout > 0.0:
                    return best_state
                
                new_st = copy.deepcopy(st)
                new_st.play(mv)
                if level <= 1:
                    new_st = self.playout(new_st, heuristic_w)
                else:
                    new_st = self.nmcs(protein, frag_mol, new_st, level - 1, heuristic_w, verbose, initial_level)
                
                combined_smiles = new_st.combinedFragSmiles(frag_mol)
                
                new_st_lip = new_st.frag_lipinskiness(combined_smiles)

                new_st_kd = new_st.frag_kdscore(protein, frag_mol)
                writeline(str(new_st.smile_to_smile(new_st.SMILE))+ " " + str(new_st_kd) +"\n", f"{self.registerName}_dock" )

                
                if (level <= initial_level -1 and new_st_kd < best_state_kd) or (level == initial_level and new_st_lip > best_state_lip): 
                    best_state = new_st
                    best_state_kd = new_st_kd
                    best_state_lip = new_st_lip
                

                    # if best_state_kd < self.best_yet:
                    #     self.best_yet = best_state_kd
                    #     # best_affinity_score = str(best_state_score-best_state.lipinskiness())
                    #     elapsed = time.perf_counter() - self.start_time
                        
                    #     writeline(str(elapsed)+ " " + str(best_state.smile_to_smile(best_state.SMILE)) +" "+ str(best_state.lipinskiness()) + " " + str(best_state_kd)+" " + "\n", f"{self.registerName}_local")
                        

            # if State.CONSIDER_NON_TERM: # early termination check
            #     if len(best_state.seq) == len(st.seq): 
            #         break
            
            if len(best_state.seq) == len(st.seq):
                break

            st.play(best_state.seq[len(st.seq)]) # st updated by playing the next move found (continues until either 'st.terminal()' or 'len(moves)==0' is met)
        
        # if State.CONSIDER_NON_TERM:
        #     return best_state
        if best_state_kd != 1000 and combined_smiles is not None:

            writeline(str(time.time() - self.start_time)+ " " + combined_smiles +  " "+ str(best_state_lip) + " " + str(best_state_kd) + "\n", f"{self.registerName}" )
            
        return st

def launch_nmcs(protein, frag_mol, init_st: State,level, heuristic_w, verbose, timeout, register_name):
    expe = frag_revNMCS()
    expe.timeout = timeout
    expe.registerName = register_name

    st = expe.nmcs(protein, frag_mol, init_st, level, heuristic_w, verbose, level)
    return st

        




    

    

        
        
