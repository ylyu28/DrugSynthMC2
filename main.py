import time
from tools import resultSaver
from models import SMILESgen
from methods import NMCS
import copy 

# def dsmc():
#     molGenState = SMILESgen.State.new()

#     if False:
#         molGenState = SMILESgen.State.make_from_string("C1=C(CCC")
    
#         print(f"Open covalence: {molGenState.nestingOpenCovalence}")
#         print(f"Open nesting ASAP: {molGenState.open_nesting_ASAP}")
#         print(f"Nesting Cycle To Close: {molGenState.nestingCycleToClose}")
#         print(f"Finish ASAP: {molGenState.finish_ASAP}")

#         start_time = time.perf_counter()

#         for m in molGenState.legal_moves():
#             print(f"atom: {m.atom}, nesting: {m.nesting}, close nesting: {m.closeNesting}, cycle: {m.cycle}, double: {m.doubleLink}")
#             print(f"prior value: {molGenState.heuristic(m).exp()}")
#             print(f"len: {len(molGenState.legal_moves())}")
#             print(f"time: {time.perf_counter() - start_time}")

#         while True:
#             pass
    
#     v = True # print reached best score when best score is reached
#     FDA_OtoC_ratio = [0.117, 0.135, 0.156, 0.139, 0.113, 0.115, 0.061, 0.037, 0.045, 0.009, 0.016, 0.013, 0.005, 0.006, 0.002, 0.003, 0.003, 0.001, 0.001, 0.013] # histogram of 20 bins, each has width = 100/20 = 5% ie ratio starts from 0.00 - 0.05 (first bin), all the way up to 0.95-1.00 (last bin, where equal count of O to C in one molecule)
#     FDA_NtoC_ratio = [0.266, 0.180, 0.184, 0.128, 0.076, 0.058, 0.031, 0.011, 0.016, 0.003, 0.016, 0.011, 0.001, 0.003, 0.000, 0.001, 0.002, 0.000, 0.000, 0.004] 
#     OtoC_bins = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
#     NtoC_bins = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]      
#     OtoC_ratio = OtoC_bins.copy()
#     NtoC_ratio = NtoC_bins.copy()

#     molecule_count = 0

#     start_time = time.perf_counter()

#     while molecule_count < 1000:

#         OtoC_diff = 0.0
#         for i in range(len(FDA_OtoC_ratio)):
#             if FDA_OtoC_ratio[i] - OtoC_ratio[i] > OtoC_diff:
#                 OtoC_diff = FDA_OtoC_ratio[i] - OtoC_ratio[i] # calculates the maximum positive difference between the target 'FDA_OtoC_ratio' to the current 'OtoC_ratio'
#                 target_OtoC_ratio = float(i) / len(FDA_OtoC_ratio) # largest diff is found for a bin with max frequency. this max frequency is set to be target ratio

#         NtoC_diff = 0.0
#         for i in range(len(FDA_NtoC_ratio)):
#             if FDA_NtoC_ratio[i] - NtoC_ratio[i] > NtoC_diff:
#                 NtoC_diff = FDA_NtoC_ratio[i] - NtoC_ratio[i]
#                 target_NtoC_ratio = float(i) / len(FDA_NtoC_ratio)
        
#         targetState = molGenState.clone()

#         st = NMCS.launch_nmcs(targetState, level=3, heuristic_w= 1.0, verbose=v, timeout=0.0, register_name="NMCS_SMILEGEN")

#         if st.reached_best_score:
#             molecule_count += 1
#             s = []
#             for c in st.SMILE:
#                 s.append(c)
            
#             if len(s) > 3 and s[-1] == ')' and s[-3] == '(':
#                 if s[-2] in ['1','2','3','4','5','6','7','8,','9']:
#                     s.remove(s[-1])
#                     s.remove(s[-2])
        
#             s2 = st.smile_to_smile(st.SMILE)
#             print(f"{s2}")
#             resultSaver.writeline(s2+"\n", "SMILES_generated/testCycles/test2")


#             carbon_count = 0.0
#             nitro_count = 0.0
#             oxygen_count = 0.0
#             for i in st.SMILE:
#                 if i == 'N':
#                     nitro_count += 1
#                 if i == 'C':
#                     carbon_count += 1
#                 if i == 'O' or i == 'U' or i == 'M':
#                     oxygen_count += 1

#             for i in range(len(OtoC_bins)):
#                 if oxygen_count/carbon_count < (i+1)/len(OtoC_bins):
#                     OtoC_bins[i] += 1
#                     break
        
#             for i in range(len(NtoC_bins)):
#                 if nitro_count/carbon_count < (i+1)/len(NtoC_bins):
#                     NtoC_bins[i] += 1
#                     break

#             sum = 0
#             for  i in range(len(OtoC_bins)):
#                 sum += OtoC_bins[i] # will give total molecule_count
#             for i in range(len(OtoC_ratio)):
#                 OtoC_ratio[i] = OtoC_bins[i]/sum
        
#             sum = 0
#             for i in range(len(NtoC_bins)):
#                 sum += NtoC_bins[i]
#             for i in range(len(NtoC_ratio)):
#                 NtoC_ratio[i] = NtoC_bins[i]/sum
        
#             v = False

#         else:
#             pass
    
#     print(f"Total time: {time.perf_counter() - start_time}")











if __name__ == '__main__':
    molGenState = SMILESgen.State.new()

    if False:
        molGenState = SMILESgen.State.make_from_string("C1=C(CCC")

        print(f"Open covalence: {molGenState.nestingOpenCovalence}")
        print(f"Open nesting ASAP: {molGenState.open_nesting_ASAP}")
        print(f"Nesting Cycle To Close: {molGenState.nestingCycleToClose}")
        print(f"Finish ASAP: {molGenState.finish_ASAP}")

        start_time = time.perf_counter()

        for m in molGenState.legal_moves():
            print(f"atom: {m.atom}, nesting: {m.nesting}, close nesting: {m.closeNesting}, cycle: {m.cycle}, double: {m.doubleLink}")
            print(f"prior value: {molGenState.heuristic(m).exp()}")
            print(f"len: {len(molGenState.legal_moves())}")
            print(f"time: {time.perf_counter() - start_time}")

        while True:
            pass

    v = True # print reached best score when best score is reached
    FDA_OtoC_ratio = [0.117, 0.135, 0.156, 0.139, 0.113, 0.115, 0.061, 0.037, 0.045, 0.009, 0.016, 0.013, 0.005, 0.006, 0.002, 0.003, 0.003, 0.001, 0.001, 0.013] # histogram of 20 bins, each has width = 100/20 = 5% ie ratio starts from 0.00 - 0.05 (first bin), all the way up to 0.95-1.00 (last bin, where equal count of O to C in one molecule)
    FDA_NtoC_ratio = [0.266, 0.180, 0.184, 0.128, 0.076, 0.058, 0.031, 0.011, 0.016, 0.003, 0.016, 0.011, 0.001, 0.003, 0.000, 0.001, 0.002, 0.000, 0.000, 0.004] 
    OtoC_bins = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    NtoC_bins = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]      
    OtoC_ratio = OtoC_bins.copy()
    NtoC_ratio = NtoC_bins.copy()

    molecule_count = 0

    start_time = time.perf_counter()

    while molecule_count < 100:

        OtoC_diff = 0.0
        for i in range(len(FDA_OtoC_ratio)):
            if FDA_OtoC_ratio[i] - OtoC_ratio[i] > OtoC_diff:
                OtoC_diff = FDA_OtoC_ratio[i] - OtoC_ratio[i] # calculates the maximum positive difference between the target 'FDA_OtoC_ratio' to the current 'OtoC_ratio'
                target_OtoC_ratio = float(i) / len(FDA_OtoC_ratio) # largest diff is found for a bin with max frequency. this max frequency is set to be target ratio

        NtoC_diff = 0.0
        for i in range(len(FDA_NtoC_ratio)):
            if FDA_NtoC_ratio[i] - NtoC_ratio[i] > NtoC_diff:
                NtoC_diff = FDA_NtoC_ratio[i] - NtoC_ratio[i]
                target_NtoC_ratio = float(i) / len(FDA_NtoC_ratio)
        
        targetState = copy.deepcopy(molGenState)

        st = NMCS.launch_nmcs(targetState, level=3, heuristic_w= 1.0, verbose=v, timeout=0.0, register_name="NMCS_SMILEGEN")

        if st.reached_best_score:
            molecule_count += 1
            s = []
            for c in st.SMILE:
                s.append(c)
            
            # if len(s) > 3 and s[-1] == ')' and s[-3] == '(':
            #     if s[-2] in ['1','2','3','4','5','6','7','8,','9']:
            #         s.remove(s[-1]) # removes ')'
            #         s.remove(s[-2]) # removes '('
        
            s2 = st.smile_to_smile(st.SMILE)
            print(f"{s2}")
            resultSaver.writeline(s2+"\n", "SMILES_generated/testCycles/test2")


            carbon_count = 0.0
            nitro_count = 0.0
            oxygen_count = 0.0
            for i in st.SMILE:
                if i == 'N':
                    nitro_count += 1
                if i == 'C':
                    carbon_count += 1
                if i == 'O' or i == 'U' or i == 'M':
                    oxygen_count += 1

            for i in range(len(OtoC_bins)):
                if oxygen_count/carbon_count < (i+1)/len(OtoC_bins):
                    OtoC_bins[i] += 1
                    break
        
            for i in range(len(NtoC_bins)):
                if nitro_count/carbon_count < (i+1)/len(NtoC_bins):
                    NtoC_bins[i] += 1
                    break

            sum = 0
            for  i in range(len(OtoC_bins)):
                sum += OtoC_bins[i] # will give total molecule_count
            for i in range(len(OtoC_ratio)):
                OtoC_ratio[i] = OtoC_bins[i]/sum
        
            sum = 0
            for i in range(len(NtoC_bins)):
                sum += NtoC_bins[i]
            for i in range(len(NtoC_ratio)):
                NtoC_ratio[i] = NtoC_bins[i]/sum
        
            v = False

        else:
            pass

    print(f"Total time: {time.perf_counter() - start_time}")
