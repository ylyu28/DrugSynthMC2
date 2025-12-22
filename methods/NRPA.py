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

    def nrpa(self, protein, rLevel, policy, initLevel):
        if rLevel == 0:
            return self.playout(State.new(),policy)
        bestLip = -1000
        bestKD = 1000
        bestState = State.new()
        for _ in range(100):
            newState = self.nrpa(protein,rLevel-1,policy,initLevel)
            newStateLip = newState.lipinskiness()
            newStateKD  = newState.kd_score(protein)
            if (rLevel <= initLevel -1 and newStateKD < bestKD) or (rLevel == initLevel and newStateLip > bestLip):
                bestLip = newStateLip
                bestKD = newStateKD
                bestState = newState
                # now if bestState gets updated, also update the policy
                # this updated policy would also be used in the nrpa for the next iteration
                policy = self.adapt(bestState.seq, policy)
        return bestState


    # def nrpa(self, protein, st:State, level, policy, initial_level):
    #     best_state = copy.deepcopy(st)
    #     best_state_lip = -1000
    #     best_state_kd = 1000
        
    #     polp = np.array(policy)

    #     while not st.terminal():
    #         moves = st.legal_moves()
    #         if len(moves)==0:
    #             break

    #         for mv in moves:

    #             new_st = copy.deepcopy(st)
    #             new_st.play(mv)

    #             for i in range(100):
    #                 if level <= 1:
    #                     new_st = self.playout(st, policy)
    #                 else:
    #                     new_st, polp = self.nrpa(protein, st, level - 1, polp, initial_level)
                    
    #                 new_st_lip = new_st.lipinskiness()
    #                 new_st_kd = new_st.kd_score(protein)
    #                 writeline(str(new_st.smile_to_smile(new_st.SMILE))+ " " + str(new_st_kd) +"\n", f"{self.registerName}_dock" )
    #                 if (level <= initial_level -1 and new_st_kd < best_state_kd) or (level == initial_level and new_st_lip > best_state_lip): 
    #                     best_state = new_st
    #                     best_state_kd = new_st_kd
    #                     best_state_lip = new_st_lip
    #                     polp = self.adapt(best_state.seq, polp)

    #         if len(best_state.seq) == len(st.seq):
    #             break

    #         st.play(best_state.seq[len(st.seq)]) 
    #         st_smile = st.smile_to_smile(st.SMILE)
    #         writeline(str(time.time() - self.start_time)+ " " + st_smile + " " + str(best_state_lip) + " "+ str(best_state_kd) + "\n", f"{self.registerName}_sc" )
        

    #     if best_state_kd != 1000 and best_state_lip != -2:
    #         smile = st.smile_to_smile(st.SMILE)
    #         writeline(str(time.time() - self.start_time)+ " " + smile + " " + str(best_state_lip) + " "+ str(best_state_kd) + "\n", f"{self.registerName}" )
            

    #     return st, polp
    
def launch_nrpa(protein, init_st: State, level, timeout, register_name):

    expe = NRPA()
    expe.timeout = timeout
    expe.registerName = register_name

    policy = [0.0] * 16 # starting policy for the seven types of legal moves

    st = expe.nrpa(protein, level, policy, level)
    return st
