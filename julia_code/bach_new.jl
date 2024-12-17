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

function enroll(values::Vector{Float64}, n::Int, m::Int, codeword::String)::Vector{Union{LinearCombination, Nothing}}
    # Assumption: input values are average-free
    std_dev = Statistics.stdm(values, 0)

    # Estimate the new distribution and guess the quantizing bounds
    dist = Normal(0, sqrt(n) * std_dev)
    quantiles = [quantile(dist, i / 2^m) for i in 1:(2^m-1)]

    weights = generate_n_bit_numbers_alpha(n, 1) 

    # Same procedure as every year 
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
            return v[1].combination[rand(filtered)]
        end

        # Return linear combination that best approximates that codeword
    end,zip(linearcombinations, par_codeword))

end
