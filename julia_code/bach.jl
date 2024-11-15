struct bach_props
    bits::Int
    bounds::Vector{Float32}
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

function greet()
    println("Hello moi")
function enroll(bounds, values::Vector{LinearCombinationSet})
    map(set -> begin
        # Filter combinations with non-negative values
        valid_combinations = filter(comb -> comb.value >= 0, set.combination)
        # Find the combination with the value closest to any bound
        max_dist_comb = findmax(comb -> minimum(abs(comb.value .- bounds)), valid_combinations)[2]
        valid_combinations[max_dist_comb]
    end, values)
end