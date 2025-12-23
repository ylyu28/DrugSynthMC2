from models.SMILESgen_cyanide import State
from tools.resultSaver import writeline
import time
import random
import math
import numpy as np
import copy



class NRPA:
    def __init__(self):
        self.timeout = -1.0
        self.registerName= ""
        self.start_time = time.perf_counter()


    def playout(self, st:State, policy):
        while not st.terminal():
            moves = st.legal_moves()
            weight = 0
            for mv in moves:
                weight = weight + math.exp(policy[st.code(mv)]+st.heuristic(mv))

            stop = random.random() * weight
            move = 0
            weight = 0
            while True:
                weight = weight + math.exp(policy[st.code(moves[move])]+st.heuristic(moves[move]))
                if weight >= stop:
                    break
                move = move + 1
            st.play(moves[move])
        return st


    def adapt(self, sequence, policy):
        """
        Policy is adapted on the best current seq found (initiating a new State to trace the best sequence)
        by increasing the weight of the best actions 
        then descreasing the weights of all the moves proportionally to their probabilities of being played
        """
        alpha = 1
        polp = np.array(policy)
        st = State.new()
        while not st.terminal():
            moves = st.legal_moves()
            z = 0
            print(st.heuristic)
            for mv in moves:
                z = z + math.exp(policy[st.code(mv)] + st.heuristic(mv))
            move = sequence[len(st.seq)] # move is the legal move chosen 
            polp[st.code(move)] += alpha

            for i in range(len(moves)):
                proba = math.exp(policy[st.code(moves[i])] + st.heuristic(moves[i]))/z
                pAdjusted = alpha * proba
                polp[st.code(moves[i])] -= pAdjusted

            st.play(move)
        
        return polp

    def nrpa(self, protein, r_level, policy, init_level):
        if r_level == 0:
            return self.playout(State.new(),policy)
        best_lip = -1000
        best_kd = 1000
        best_state = State.new()
        for _ in range(100):
            new_state = self.nrpa(protein,r_level-1,policy,init_level)
            new_state_lip = new_state.lipinskiness()
            new_state_kd  = new_state.kd_score(protein)
            writeline(str(new_state.smile_to_smile(new_state.SMILE))+ " " + str(new_state_kd) +"\n", f"{self.registerName}_dock" )
            if (r_level <= init_level -1 and new_state_kd < best_kd) or (r_level == init_level and new_state_lip > best_lip):
                best_state_lip = new_state_lip
                best_state_kd = new_state_kd
                best_state = new_state
                if best_state_kd != 1000 and best_state_lip != -2:
                    smile = new_state.smile_to_smile(new_state.SMILE)
                    writeline(str(time.time() - self.start_time)+ " " + smile + " " + str(best_state_lip) + " "+ str(best_state_kd) + "\n", f"{self.registerName}" )

                # now if bestState gets updated, also update the policy
                # this updated policy would also be used in the nrpa for the next iteration
                policy = self.adapt(best_state.seq, policy)

        return best_state


    
def launch_nrpa(protein, level, timeout, register_name):

    expe = NRPA()
    expe.timeout = timeout
    expe.registerName = register_name

    policy = [0.0] * 16 # starting policy for the seven types of legal moves

    st = expe.nrpa(protein, level, policy, level)
    return st
