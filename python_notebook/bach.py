# For every entry in lincomb_sets, iterate over every possible linear combination and calulate the distance to the bounds. 
# Find the nearest bound and save the bound for that linear combination. 
# Then, select the linear combinations with the highest distance, the resulting
# weights, is the thing that we want to have  

# Use parallel computing to speed up the process 
# lincomb_set has every possible value of linear combinations for a set of input values

from functools import partial

def enroll(lincomb_sets, bounds):
    # Iterate over every linear combination set 
    # TODO let this run in parallel to increase performance
    # TODO Way to inefficient, needs to be optimized here
    optimal_combinations = []

    for set in lincomb_sets: 
        possible_optimal_combinations = []
        for combination in set.combination:
            # Distances to every bound for a linear combination
            distances = []
            for bound in bounds: 
                distances.append(abs(combination.value - bound))
            possible_optimal_combinations.append((combination, min(distances)))
        #print(possible_optimal_combinations[0])
        optimal_combinations.append(max(possible_optimal_combinations, key=lambda x: x[1]))
        
    return optimal_combinations


def in_range(combination, nogozone):
    return nogozone[0] <= combination.value <= nogozone[1]


# Didnt really work..
def enroll_unfavour(lincomb_sets, bounds, width):
    # Generate an area arount the 0 that is considered unfavourable and take a look at what happens
    
    # Define the range
    nogozone = [-width/2, width/2]

    optimal_combinations = []

    for set in lincomb_sets:
        possible_optimal_conbinations = []
        for combination in set.combination: 
            distances = []
            for bound in bounds:
                distances.append(abs(combination.value - bound))
            possible_optimal_conbinations.append((combination, min(distances)))
        # Everything until here stays the same. Now we will check if in the possible_optimal_combinations for a value inside the 
        # range
        in_range_partial = partial(in_range, nogozone=nogozone)
        possible_optimal_combinations_filtered = filter(
            lambda x: not in_range_partial(x[0]),
            possible_optimal_conbinations
        )
        optimal_combinations.append(max(possible_optimal_combinations_filtered, key=lambda x: x[1]))
    return optimal_combinations