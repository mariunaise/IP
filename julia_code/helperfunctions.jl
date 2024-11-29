"Returns a list of n-bit numbers where 0 is exchanged with -a and 1 with a"
function generate_n_bit_numbers_alpha(n, a)
    numbers = []
    for i in 0:(2^n - 1)
        binary_str = bitstring(i)[end-n+1:end]
        transformed = [c == '0' ? a : -a for c in binary_str]
        push!(numbers, transformed)
    end
    return numbers
end

function create_linearcombinations(inputs, weights, n) 
    collect(
    map(
        set -> begin
            LinearCombinationSet(
                collect(map(
                    weights -> begin
                        LinearCombination(weights, set, sum(weights .* set))
                    end,
                    weights
                ))
            )
        end, 
        Iterators.partition(inputs, n)
    )
)
end

function find_interval(A::Float32, B::Vector{Float64})
"""
This function returns the two entries of the vector B such that A is between them. 
If A is smaller than the smallest entry, lower is nothing and vice versa.
"""
    # Sort the bounds vector in ascending order
    sorted_B = sort(B)
    # Initialize the two return values
    lower = nothing
    upper = nothing 

    for i in eachindex(sorted_B)
        if A < sorted_B[i]
            upper = sorted_B[i]
            lower = i > 1 ? sorted_B[i-1] : nothing
            break
        elseif A == sorted_B[i]
            lower = sorted_B[i]
            upper = sorted_B[i]
            break
        end

        if A > sorted_B[end]
            lower = sorted_B[end]
        end
    end

    # Only return two elements if both of them contain a value, otherwise only return one

    if lower === nothing
        return upper
    elseif upper === nothing
        return lower
    end
    
    return (lower, upper)
end