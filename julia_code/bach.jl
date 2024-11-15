include("helperfunctions.jl")

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

# Filter Bach Input only naively based on the maximum distance of some values 
function bach_filter_naive(props::bach_props, values::Vector{LinearCombinationSet})
    # First Iteration: Optimize the values away from the 0 to obtain the 
    # filterd values for 1 bit quanitzation 
    

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

# Starter function, optimize values away form the mean of the original distribution
function start(values::Vector{LinearCombinationSet})
    map(set -> begin    
        set.combination[findmax(comb -> abs(comb.value), set.combination[2])]
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

# Function takes the set of bounds generated in the step before and filteres away 
function splitter(bounds, values::Vector{LinearCombinationSet})

end

function bach(inputs::Vector{Float64})
    # Initialize with mean of the distribution as first bound
    bmv = [mean(inputs)]
    print(bmv)
end