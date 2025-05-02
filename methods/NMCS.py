from models.SMILESgen import State
from tools.calc import softmaxChoice
from tools.resultSaver import writeline
import time
import random


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
        best_state = st.clone()
        best_state_score = 0.0

        if State.CONSIDER_NON_TERM or st.terminal():
            best_state_score = best_state.score()

        while not st.terminal():
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
                sc = st.score()
                if sc > best_state_score:
                    best_state_score = sc
                    best_state = st.clone()

        return best_state if State.CONSIDER_NON_TERM else st
    

    def nmcs(self, st: State, level, heuristic_w, verbose: bool):
        best_state = st.clone()
        best_state_score = -1.0

        while not st.terminal(): # runs until the state is terminal or no legal moves are left
            moves = st.legal_moves()
            if len(moves) == 0:
                break
            for mv in moves:
                if (time.time() - self.start_time) > self.timeout and self.timeout > 0.0:
                    return best_state
                
                new_st = st.clone()
                new_st.play(mv)
                if level <= 1:
                    new_st = self.playout(new_st, heuristic_w)
                else:
                    new_st = self.nmcs(new_st, level - 1, heuristic_w, verbose)
                new_st_score = new_st.score()

                if new_st.reached_best_score:
                    if verbose: 
                        print(f"Reached best score")
                    return new_st
                
                if new_st_score > best_state_score: # update best state if new state is better
                    best_state = new_st
                    best_state_score = new_st_score
                    if best_state_score > self.best_yet:
                        self.best_yet = best_state_score
                        elapsed = time.perf_counter() - self.start_time
                        writeline(str(elapsed)+ " " + str(best_state_score) + "\n", self.registerName.clone())

            if State.CONSIDER_NON_TERM: # early termination check
                if len(best_state.seq) == len(st.seq): # if score not improved in this iteration, consider termination at non_terminal state
                    break
            
            st.play(best_state.seq[len(st.seq)]) # st updated by playing the next move found (continues until either 'st.terminal()' or 'len(moves)==0' is met)
        
        if State.CONSIDER_NON_TERM:
            return best_state
        
        return st

def launch_nmcs(init_st: State, level, heuristic_w, verbose, timeout, register_name):
    expe = NMCS()
    expe.timeout = timeout
    expe.registerName = register_name

    st = expe.nmcs(init_st, level, heuristic_w, verbose)
    return st 

        




    

    

        
        
