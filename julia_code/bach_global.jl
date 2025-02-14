include("helperfunctions.jl")

using Statistics
using Distributions

mutable struct LinearCombination
    """
    Struct to represent a linear combination of input values. 
    Contains a weights field with the weights used for the linear combination, a helperdatabits field that contains the information which input values are added and which are subtracted, an inputs field that contains the input values used for the linear combination, a value field that contains the result of the linear combination.
    """
    weights::Vector{Float32}
    # Per convention of julia, a false value means addition, a true one subtraction
    helperdatabits::Vector{Vector{Bool}}
    inputs::Vector{Float32}
    value::Float32
end 

struct LinearCombinationSet
    """
    Struct to represent a set of linear combinations. 
    Used by the enrollment function to store multiple linear combinations that are possible solutions for one set of input values.
    """
    combination::Vector{LinearCombination}
end 

function enroll_optimized(linearcombinations, bounds) # 340ms
    """
    Enrollment function 
    Takes a Vector of LinearCombinationSet objects as possible input combinations and a vector of bounds to perform the sieving operation. 
    Bounds has to be a vector of float values that are symmetrically distributed around 0.
    """
    # Out of every LinearCombinationSet find the LinearCombination that best optimizes away from each bound 
    result = map(set -> begin 
        # Return the LinearCombination object that is furthest away from any of the bounds
        distances = map(combination -> begin 
            map(bound -> begin 
                abs(combination.value - bound) 
            end, bounds)  
        end, set.combination)  
            
        return set.combination[argmax(map(v -> minimum(v), distances))] 
            
    end,linearcombinations)
    
    return map(comb -> comb.value, result)

end

function enroll(values, n, m, iterations, alpha, distfactor, boundsparam=nothing)::Tuple{Vector{LinearCombination}, Vector{Float32}}
    """
    DEPRECATED FUNCTION, USE enroll_optimized INSTEAD
    Function to perform a global enrollment. 
    alpha is there to limit the distance a bound can travel in one iteration. 
    distfactor is a scalar factor to estimate the new distribution.
    boundsparam is an optional parameter to set the bounds manually.
    iterations is the number of iterations the algorithm should run.
    Returns the LinearCombinations that are the solution of the last iteration.
    """
    # Generate Weights 
    weights = generate_n_bit_numbers_alpha(n, 1)

    # Create LinearCombinationSets 
    linearcombinations = create_linearcombinations(values, weights, n)

    # Initial bounds are quantiles of the n-fold convolution of input values
    #std_dev = Statistics.stdm(map(set -> begin 
    #    map(comb -> comb.value, set.combination)
    #end, linearcombinations), 0.0)
    
    std_dev = 1

    # Estimate the new distribution and guess the quantizing bounds
    dist = Normal(0, distfactor * std_dev)

    if boundsparam == nothing
        bounds = [quantile(dist, i / 2^m) for i in 1:(2^m-1)]
    else
        bounds = boundsparam
    end
    
    println("Old bounds: ", bounds)

    solution_combinations = 0

    # Iterate multiple times and update the linearcombinations that exist
    for i in range(1, iterations)
        # Out of every LinearCombinationSet find the LinearCombination that best optimizes away from each bound 
        part_solution = map(set -> begin 
            # Return the LinearCombination object that is furthest away from any of the bounds
            distances = map(combination -> begin 
                map(bound -> begin 
                    abs(combination.value - bound) 
                end, bounds)  
            end, set.combination)  
            
            return set.combination[argmax(map(v -> minimum(v), distances))] 
            
        end,linearcombinations)
        
        # Estimate new quantizing bounds based on the resulting LinearCombination values. 
        # Limit the moving of the bounds using the mean of the x value of the bound before and the one after the optimization step 
        sorted_values = sort(map(comb -> comb.value, part_solution))
        new_bounds = ecdf_bounds(sorted_values, m) 

        # Calculate the distance each bound travelled to its new location. 
        # Limit that distance with a scalar factor
        bound_distances = map(bnds -> bnds[1] - bnds[2] , zip(bounds, new_bounds))
        println("The distance to bounds in iteration", i ,"are :", bound_distances)
        # First bound travelled bound_distances[1] into a direction. 
        # Multiply the distance with a scalar factor alpha to limit the distance a bound travells 

        #println(typeof(bound_distances))

        bound_distances = bound_distances .* alpha

        # Update the bounds with new bound distances 
        updated_bounds = [bounds[i] + bound_distances[i] for i in 1:length(bounds)]
        println("New bounds: ", updated_bounds)

        bounds = updated_bounds

        solution_combinations = part_solution
    end
    
    # At the end, return the part_solution of the last iteration 
    return (solution_combinations, bounds)
end

function reconstruct(values::Vector{Float64}, n, weights::Vector{Vector{Float32}})
    """
    Function to reconstruct the LinearCombination values from found weights 
    Needs the (normal distributed) input values, the number of addends n and the vector of weights for every linear combination as two dimensional array.
    """
    map(part -> sum(part[1] .* part[2]), zip(Iterators.partition(values,n), weights)) 
end
