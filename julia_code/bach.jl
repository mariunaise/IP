include("helperfunctions.jl") 

ENV["BROWSER"] = "loupe"
using Gadfly
using Statistics
using Distributed
using Base.Iterators: flatten

#struct bach_props
#    bits::Int
#    bounds::Vector{Float32}
#end

mutable struct LinearCombination
    weights::Vector{Float32}
    # Per convention of julia, a false value means addition, a true one subtraction
    helperdatabits::Vector{Vector{Bool}}
    inputs::Vector{Float32}
    value::Float32
end 

struct LinearCombinationSet
    combination::Vector{LinearCombination}
end 

function bach_enroll_first(values::Vector{LinearCombinationSet})
    # Iterate through every set of LinearCombinationSet and only choose the LinearCombination 
    # for which the absolute value of the value field is the biggest
    map(set -> begin
        # findmax reutrns a tuple of the maximum value and its index in the collection
        # We want to access the index of the maximum based on the lambda function
        # that computes the absolute value of a given linear combiantion result 
        set.combination[findmax(comb -> abs(comb.value), set.combination)[2]]
    end, values) 
end

function enroll(bounds, values::Vector{LinearCombinationSet})
    map(set -> begin
        # Filter combinations with non-negative values
        valid_combinations = filter(comb -> comb.value >= 0, set.combination)
        # Find the combination with the value closest to any bound
        max_dist_comb = findmax(comb -> minimum(abs(comb.value .- bounds)), valid_combinations)[2]
        valid_combinations[max_dist_comb]
    end, values)
end

function enroll2(bounds, values::Vector{LinearCombinationSet})
    map(set -> begin
    #println("TYPE OF VALID COMBINATIONS" * string(typeof(set.combination)))
        valid_combinations = set.combination
        max_dist_comb = findmax(comb -> minimum(abs(comb.value .- bounds)), valid_combinations)[2]
        valid_combinations[max_dist_comb]
    end, values)
end

# Starter function, optimize values away form the mean of the original distribution
function start(values::Vector{LinearCombinationSet})
    """
    Starts the optimization process by maximizing the distance to the point 0
    Returns a symmetrical distribution containing two skewed normal distributions
    """
    map(set -> begin    
        set.combination[findmax(comb -> abs(comb.value), set.combination)[2]]
    end,values)
end

function filterer(bounds, initialCombination::LinearCombination, tryingCombinations::LinearCombinationSet)::LinearCombinationSet
    # Iterate over every possible new linear combination and filter the ones that step over the bounds
    valid_combinations =  filter(trialCombination -> begin
        all(bound -> signbit(initialCombination.value - bound) == signbit(trialCombination.value - bound), bounds)
    end,tryingCombinations.combination)

    # Return a LinearCombinationSet object that contains only the valid LinearCombination objects 
    LinearCombinationSet(collect(map(combination -> begin
        combination
    end,valid_combinations)))
end

function optimizer(values::Vector{LinearCombination}, bounds, fraction, n)
    # First, define the fractional weight vector to be added to the original weights 
    frac_weights = generate_n_bit_numbers_alpha(n, fraction)

    # Find median of the input values to define the next bound value
    median_val = Statistics.median(map(x -> x.value, values))

    # Create new Linear Combination sets with the fractional weight
    # This variable is a Vector::LinearCombinationSet again
    new_combinations = collect(map(combination -> begin
    LinearCombinationSet(collect(map(weights_vector -> begin
        weights = combination.weights + weights_vector
        helperdata = map(v -> vcat(v[1], signbit(v[2])), zip(combination.helperdatabits, weights_vector))
            return LinearCombination(weights, helperdata ,combination.inputs, sum(weights .* combination.inputs))
        end,frac_weights)))
    end,values))

    # zip together the initial combinations and the LinearCombinationSets so the resutling structure 
    # is a Vector of tuples Vector::(intiialCombination, trialCombinationSet)
    zipped_combinations = collect(zip(values, new_combinations))

    # Now somehow magically filter all the values that cross a quantizing bound 
    
    # 1. Find out which quantizing bounds are relevant for this specific entry 
    # This can be achieved by taking a look at the median value and finding the interval
    # for which a is eigher bigger, smaller on inbetween one or two entries of the bounds vector 
    relevant_bounds = find_interval(median_val, bounds)

    # 2. With the bounds found, filter the Values of the linear combination sets that exceed a quantizing bound 
    # after the fractional weight addition 

    # To achieve this, we will call a filter function with the bounds, the initial combination before the addition of the 
    # fractional weight and the set of new weights to filter the set of new weights accordingly
    filtered_combinations = collect(map(zipped_combination -> begin 
        filterer(relevant_bounds, zipped_combination[1], zipped_combination[2])
    end,zipped_combinations))

    # 3. After the values have been filtered, perform the usual boundary distance maximization stategy as always to 
    # get the splittted Gumbel distribtuion of this set.

    # Do the magic here 
    result = enroll2(median_val, filtered_combinations)

    # Add the new median value to the bounds vector
    push!(bounds, median_val)

    return result
    
end

function splitter(bounds, values::Vector{LinearCombination})
    result = []
    # Iterate through every interaval of the bounds vector and create a new distribution with each interval
    for bound in bounds 
       push!(result, filter(comb -> comb.value <= bound, values))
       filter!(comb -> comb.value > bound, values)
    end
    push!(result, values)
    return result
end

# BACH function to shape the input values for m bit quantization using a linear combination of n addends 
function bach(inputs::Vector{Float64}, n, m, number_sequence)::Tuple{Vector{LinearCombination}, Vector{Float32}}
    """
    # BREAK (Boundary Recursive Epic Adaptive (Adventure) Klustering)
    """
    # Initial weights are 1 and -1
    weights = generate_n_bit_numbers_alpha(n, 1)

    # Create initial linear combinations and optimize them away from the 0 
    # linear_combinations here contians already optimized values and weights for the 1 bit quanization and can be used here if m = 1
    linear_combinations = start(create_linearcombinations(inputs, weights, n))

    # Define a vector with quantizin bounds and initialize it with 0
    quantizing_bounds = [0.0]

    # During every iteration, spawn 2^(m-1) child processes optimizing their respective sub distribution with a fractional weight
    for i in 2:m 
        # Define the sub distributions based on the bounds vector
        sub_distributions = splitter(quantizing_bounds, linear_combinations)
        linear_combinations = Vector{LinearCombination}()
        shifting_index = number_sequence(i)
        for distribution in sub_distributions
            # Sort the quantizing_bounds vector in ascending order: 
            sort!(quantizing_bounds)

            # Plot the subdistributions 
            i = findfirst(item -> item == distribution, sub_distributions)
            
            # Apply the optimizing function to each sub distribution and update linear_combinations for the next iteration
            linear_combination = optimizer(distribution, quantizing_bounds, shifting_index, n)
            append!(linear_combinations, linear_combination)
            
        end
    end
    return (linear_combinations, quantizing_bounds)
end

function reconstruct_weights(helperdata, n, m, number_sequence)
    # Needs to return a vector of weights that can be zipped together 
    # We can deduce the number sequence based on the amounts of bits we want to quantize. 
    # For every bit we want to quantize, we will need one more iteration -1 of the number_sequence function 

    # The first weight is always 1, so we will need m-1 iteration of the number_sequence function
    fractional_weights = append!([1.0], map(v -> number_sequence(v),collect(range(2, m))))

    # Based on number of addends and number of quantized bits, the helperdata vector consists of m x n matrices. 
    # TODO Montag
end

function reconstruct(inputs::Vector{Float64}, n, m, number_sequence)
      
end
