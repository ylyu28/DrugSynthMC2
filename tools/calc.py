import numpy as np
import random
import math

# Implement softmax, Manhattan distance, distance matrix

def softmaxChoice(l):
    r = random.random()
    total = sum(math.exp(x) for x in l)
    
    cumulative_prob = 0
    for i, value in enumerate(l):
        cumulative_prob += math.exp(value) / total
        
        if cumulative_prob >= r:
            return i
        
    # if no i was returned, move is the last 
    print(f" {r=} {total=} {cumulative_prob=}")
    print(f"List values: {l}")
    # return len(l)
    return len(l) - 1 # i think we should use len(l)-1 instead

    