include("helperfunctions.jl")

using Distributions
using Statistics

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

function enroll(values::Vector{Float64}, n::Int, m::Int, codeword::String, uniform::Bool)::Vector{Union{LinearCombination, Nothing}}
    # Assumption: input values are average-free
    std_dev = Statistics.stdm(values, 0)

    # Estimate the new distribution and guess the quantizing bounds
    dist = Normal(0, std_dev)
    println("Calculating " * string(2^m) * " quantiles")

    if uniform == true
        quantiles = [i / 2^m for i in 1:(2^m)]
    else
        quantiles = [quantile(dist, i / 2^m) for i in 1:(2^m-1)]
    end

    #quantiles = [quantile(dist, i / 2^m) for i in 1:(2^m-1)]
    
    println(quantiles)

    weights = generate_n_bit_numbers_alpha(n, 1) 

    # We can also add some more weight combinations to this list to maybe increase accuracy ..
    # 
    additional_weights = 

    # Same procedure as every year 
    #linearcombinations = create_linearcombinations(values, weights, n) 

    linearcombinations = create_linearcombinations(values, weights, n)

    # Partition the codeword based on the number of bits we want to extract and convert them to integers
    par_codeword = collect(map(v -> parse(Int, v, base=2),Iterators.partition(codeword, m)))

    return map( v -> begin 
        
        # Quantize every linear combination based on the quantiles vector and return the one that 
        # (1) lies in the intended quantizing bin
        # (2) is the furthest away from any quantizing bound  
        
        # Extract values from the LinearCombination objects (will still be in the same order)
        values = map(w -> w.value ,v[1].combination)

        # Quantize the values based on the quantiles found for the input distribution
        quantized = collect(searchsortedlast.(Ref(quantiles), values))

        # Determine the indices of LinearCombinations that quantize to the par_codeword for this iteration 
        filtered = findall(x -> x == v[2] ,quantized)
   
        if isempty(filtered)
            return
        else
            #println("Thingy is not empty with " * string(length(filtered)) * "elements")
            # Determine the quantizing bounds of one of the elements in the filtered list to find an optimal solution
            relevant_bounds = find_interval(v[1].combination[1].value, quantiles)

            # Calculate the distance to every bound of the v[1].combination vector 
            # Should be a Vector{Vector{Float64}}
            distances = map(combination -> begin 
                map(bound -> begin 
                    abs(combination.value - bound) 
                end,relevant_bounds)  
            end, v[1].combination)

            return v[1].combination[argmax(map(v -> minimum(v), distances))]
            #return v[1].combination[rand(filtered)]
        end

        # Return linear combination that best approximates that codeword
    end,zip(linearcombinations, par_codeword))

end
