include("helperfunctions.jl") 

ENV["GDK_BACKEND"] = "x11"
using Gadfly

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
    println("ICH BIN HIER!!!")
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

# Starter function, optimize values away form the mean of the original distribution
function start(values::Vector{LinearCombinationSet})
    map(set -> begin    
        set.combination[findmax(comb -> abs(comb.value), set.combination)[2]]
    end,values)
end

"### Optimizing Function
Based on the amount of bounds passed, the initial filtering differs
1. If bounds only contains one value, exclude every resulting combination 
2. If bounds contains two values, sort them in ascending order and filter every resulting combination outside of these bounds.
" 
function optimizer(values::Vector{LinearCombinationSet}, bounds, fraction, n)
    # First, define the fractional weight vector to be added to the original weights 
    frac_weights = generate_n_bit_numbers_alpha(n, fraction)

end

#  
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
function bach(inputs::Vector{Float64}, n, m)
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
        for distribution in sub_distributions
            # Plot the subdistributions 
            i = findfirst(item -> item == distribution, sub_distributions)
            display(plot(x=collect(map(comb -> comb.value, distribution)), Geom.histogram(bincount=1000), Guide.title("Sub Distribution " * string(i))))
            
            # Optimize the sub distributions 
        end
    end

end