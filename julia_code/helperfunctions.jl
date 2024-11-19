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