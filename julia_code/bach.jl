include("helperfunctions.jl") 

ENV["BROWSER"] = "loupe"
using Gadfly
using Statistics
using Distributed
using Base.Iterators: flatten

struct bach_props
    bits::Int
    bounds::Vector{Float32}
end

mutable struct LinearCombination
    weights::Vector{Float32}
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
            return LinearCombination(weights, combination.inputs, sum(weights .* combination.inputs))
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


function find_bounds()
"""
Function to find the bounds 

Takes some kind of input distribution and the amout of bits to be quanitzed and returns a vector of quantizing bounds 
so the enrollment function does not have to caluclate the median of the sub distribtuions anymore but only a scalar factor 
is to be found
"""
end

# BACH function to shape the input values for m bit quantization using a linear combination of n addends 
function bach(inputs::Vector{Float64}, n, m)
    """
    # BREAK (Boundary Recursive Epic Adaptive (Adventure) Klustering)
    """
    # First, we will optimize away from the 0, so the start() function is called for that 

    #display(plot(x=inputs, Geom.histogram(bincount=1000)))

    # Initial weights are 1 and -1
    weights = generate_n_bit_numbers_alpha(n, 1)

    # Create initial linear combinations and optimize them away from the 0 
    # linear_combinations here contians already optimized values and weights for the 1 bit quanization and can be used here if m = 1
    linear_combinations = start(create_linearcombinations(inputs, weights, n))

    # Plot the result of the first iteration here
    display(plot(x=collect(map(comb -> comb.value, linear_combinations)), Geom.histogram(bincount=1000), Guide.title("Iteration 1")))

    # Define a vector with quantizin bounds and initialize it with 0
    quantizing_bounds = [0.0]

    # During every iteration, spawn 2^(m-1) child processes optimizing their respective sub distribution with a fractional weight
    for i in 2:m 
        # Define the sub distributions based on the bounds vector
        sub_distributions = splitter(quantizing_bounds, linear_combinations)
        linear_combinations = Vector{LinearCombination}()
        for distribution in sub_distributions
            # Sort the quantizing_bounds vector in ascending order: 
            sort!(quantizing_bounds)

            # Plot the subdistributions 
            i = findfirst(item -> item == distribution, sub_distributions)
            display(plot(x=collect(map(comb -> comb.value, distribution)), Geom.histogram(bincount=1000), Guide.title("Sub Distribution " * string(i))))
            
            # Apply the optimizing function to each sub distribution and update linear_combinations for the next iteration
            linear_combination = optimizer(distribution, quantizing_bounds, 0.3, n)
            display(plot(x=collect(map(comb -> comb.value, linear_combination)), Geom.histogram(bincount=1000), Guide.title("Result after optimizing sub distribution " * string(i))))

            println("bounds vector is currently: " * string(quantizing_bounds))

            append!(linear_combinations, linear_combination)
            
        end
    end

    display(plot(x=collect(map(comb -> comb.value, linear_combinations)), Geom.histogram(bincount=1000), Guide.title("Finished Result")))
end