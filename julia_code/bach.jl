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
        set.combination[findmax(comb -> abs(comb.value), set.combination)[2]]
    end, values) 
end

function greet()
    println("Hello moi")
end