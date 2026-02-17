from models.SMILESgen_cyanide import State
from tools.calc import softmaxChoice
from tools.resultSaver import writeline
import time
import random
import copy


class lazyNMCS:
    def __init__(self):
        self.best_yet = -1000.0
        self.timeout = -1.0
        self.registerName= ""
        self.start_time = time.perf_counter()
    

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
    

    def lnmcs(self, protein, st: State, level, heuristic_w, verbose: bool, initial_level, r_lip, mean_playout_count):
        """

        Tunable parameters: 
        initial_level: nested levels
        r_lip: ratio for minimum drug-likeness threshold
        mean_playout_count: number of simulated playouts for score evaluation

        """
        best_state = copy.deepcopy(st)
        best_state_lip = -1000
        best_state_kd = 1000


        while not st.terminal(): # runs until the state is terminal (no legal moves are left)
            moves = st.legal_moves()
            if len(moves) == 0:
                break
            lip_scores = [0 for m in moves]

            for i in range (len(moves)):
                mv = moves [i]

                for j in range (mean_playout_count):
                    new_st = copy.deepcopy(st)
                    new_st.play(mv)
                    new_st = self.playout(new_st, heuristic_w)
                    new_st_lip = new_st.lipinskiness()
                    lip_scores[i] += new_st_lip

            inferior = min(lip_scores)

            
            for i in range (len(moves)):
                mv = moves[i]
                if lip_scores[i] >= r_lip * inferior: # prune nodes that do not exceed the threshold
                    if (time.time() - self.start_time) > self.timeout and self.timeout > 0.0:
                        return best_state
                
                    new_st = copy.deepcopy(st)
                    new_st.play(mv)
                    if level <= 1:
                        new_st = self.playout(new_st, heuristic_w)
                    else:
                        new_st = self.lnmcs(protein, new_st, level - 1, heuristic_w, verbose, initial_level, r_lip, mean_playout_count)

                    new_st_lip = new_st.lipinskiness()

                    new_st_kd = new_st.kd_score(protein)
                    writeline(str(new_st.smile_to_smile(new_st.SMILE))+ " " + str(new_st_kd) +"\n", f"{self.registerName}_dock" )

                
                    if (level <= initial_level -1 and new_st_kd < best_state_kd) or (level == initial_level and new_st_lip > best_state_lip): 
                        best_state = new_st
                        best_state_kd = new_st_kd
                        best_state_lip = new_st_lip
                
            if len(best_state.seq) == len(st.seq):
                break

            st.play(best_state.seq[len(st.seq)]) # st updated by playing the next move found (continues until either 'st.terminal()' or 'len(moves)==0' is met)
            st_smile = st.smile_to_smile(st.SMILE)
            writeline(str(time.time() - self.start_time)+ " " + st_smile + " " + str(best_state_lip) + " "+ str(best_state_kd) + "\n", f"{self.registerName}_sc" )
        
        # if State.CONSIDER_NON_TERM:
        #     return best_state
        if best_state_kd != 1000 and best_state_lip != -2:
            smile = st.smile_to_smile(st.SMILE)
            writeline(str(time.time() - self.start_time)+ " " + smile + " " + str(best_state_lip) + " "+ str(best_state_kd) + "\n", f"{self.registerName}" )
            

        return st


def launch_lnmcs(protein, init_st: State,level, heuristic_w, verbose,timeout, r_lip, mean_playout_count, register_name):
    expe = lazyNMCS()
    expe.timeout = timeout
    expe.registerName = register_name

    st = expe.lnmcs(protein, init_st, level, heuristic_w, verbose, level, r_lip,mean_playout_count)
    return st
