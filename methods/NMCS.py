from models.SMILESgen import State
from tools.calc import softmaxChoice
from tools.resultSaver import writeline
import time
import random
import copy
import math


class NMCS:
    def __init__(self):
        self.best_yet = 0.0
        self.timeout = -1.0
        self.registerName= ""
        self.start_time = time.perf_counter()
    
    @classmethod
    def new(cls):
        return cls()
    

    def playout(self, st: State, heuristic_w): # playout, which employs softmax, used in level=1 nmcs
        best_state = copy.deepcopy(st)
        best_state_score = 0.0
        best_state_kd = 1000

        if State.CONSIDER_NON_TERM or st.terminal():
            best_state_kd, best_state_score = best_state.score()

        while not st.terminal(): # keep adding 'artificial moves'
            moves = st.legal_moves()
            if len(moves) == 0:
                break
            
            i = int(len(moves) * random.random())
            if heuristic_w != 0.0:
                weights = [heuristic_w * st.heuristic(m) for m in moves] # proability of each move
                i = softmaxChoice(weights) # choosing the move the softmax way

            mv = moves[i]
            st.play(mv) 

            if State.CONSIDER_NON_TERM:

                kd, sc = st.score()
                if sc > best_state_score:
                    best_state_score = sc
                    best_state = copy.deepcopy(st) 

        return best_state if State.CONSIDER_NON_TERM else st 
    

    def nmcs(self, st: State, level, heuristic_w, verbose: bool):
        best_state = copy.deepcopy(st)
        best_state_score = -1.0
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
                    new_st = self.nmcs(new_st, level - 1, heuristic_w, verbose)
                new_st_kd, new_st_score = new_st.score()
                writeline(str(new_st.smile_to_smile(new_st.SMILE))+ " " + str(new_st_kd) +"\n", f"{self.registerName}_dock" )

                if new_st.reached_best_score:
                    if verbose: 
                        print(f"Reached best score")
                        st_smile = st.smile_to_smile(st.SMILE)
                        # new_st_kd = new_st.kd_score
                        elapsed = time.perf_counter() - self.start_time
                        writeline(str(time.time() - self.start_time)+ " " + st_smile + " " + str(new_st_score) + " "+ str(new_st_kd) + "\n", f"{self.registerName}_sc" )
                    return new_st
                
                if new_st_score > best_state_score: # update best state if new state is better
                    best_state = new_st
                    best_state_kd = new_st_kd
                    best_state_score = new_st_score

                    if best_state_score > self.best_yet:
                        self.best_yet = best_state_score
                        # best_affinity_score = str(best_state_score-best_state.lipinskiness())
                        elapsed = time.perf_counter() - self.start_time
                        print(best_state_score)
                        writeline(str(elapsed)+ " " + str(best_state.smile_to_smile(best_state.SMILE)) +" "+ str(best_state_score) + " " + str(best_state_kd)+" " + "\n", f"{self.registerName}_local")
                        

            if State.CONSIDER_NON_TERM: # early termination check
                if len(best_state.seq) == len(st.seq): # if score not improved in this iteration, consider termination at non_terminal state
                    break
            
            st.play(best_state.seq[len(st.seq)]) # st updated by playing the next move found (continues until either 'st.terminal()' or 'len(moves)==0' is met)
            st_smile = st.smile_to_smile(st.SMILE)
            writeline(str(time.time() - self.start_time)+ " " + st_smile + " " + str(best_state_score) + " "+ str(best_state_kd) + "\n", f"{self.registerName}_sc" )
        
        if State.CONSIDER_NON_TERM:
            return best_state
        
        # writeline(str(time.time()-self.start_time)+ " " + str(st.SMILE)+ " " + str(st.score()) +"\n", "scoreMonitor" )
        
        return st

def launch_nmcs(init_st: State, level, heuristic_w, verbose, timeout, register_name):
    expe = NMCS()
    expe.timeout = timeout
    expe.registerName = register_name

    st = expe.nmcs(init_st, level, heuristic_w, verbose)
    return st 

        




    

    

        
        
