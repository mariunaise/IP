using Pkg
Pkg.activate("env")

addprocs(12)

using Distributions
using Random
@everywhere using HypothesisTests
@everywhere using CSV, DataFrames
@everywhere using Base.Threads
@everywhere using StatsBase

@everywhere include("helperfunctions.jl")

# General setup of input values
println("--Welcome to the best Brute Forcing script in the wild west!--")
println("Using ", Threads.nthreads(), " threads")

dist_mean = 0
std_dev = 1

n = 3
m = 3

dist = Normal(dist_mean, std_dev)
data = rand(dist, 1000000)

# Define bound combinations and construct array with all possible bounds that we want to test

A = collect(0.1:0.01:1)
B = collect(1.1:0.01:2)
C = collect(2.1:0.01:3)

bound_combinations = collect(Iterators.product(A, B, C))

possible_bounds = reshape(map(tuple -> begin return [-tuple[3], -tuple[2], -tuple[1], 0.0, tuple[1], tuple[2], tuple[3]] end,bound_combinations), (:))

#possible_bounds = map(value -> [-value, 0.0, value], A)

# Define lock against parallel processing race conditions

csv_lock = ReentrantLock()

# Define weights 

bach_weights = generate_n_bit_numbers_msb(n)

# Define the possible linearcombinations now 
linearcombinations = create_linearcombinations(data, bach_weights, n)

# Define a function to append the result to the CSV
@everywhere function save_to_csv(data)
    #lock(csv_lock) do
        df = DataFrame(reshape(data, 1, :), :auto)        # TODO
        CSV.write("simulation_results_3bit.csv", df, append=true)
    #end
end

@everywhere function chi_squared_test(input_values, bounds) # 6ms
    # Quantize the values based on given bounds
    quantized_values = searchsortedlast.(Ref(bounds), input_values)
    observed_frequencies = collect(DataFrames.values(countmap(quantized_values)))
    #println(observed_frequencies)
    total_values = length(input_values)
    num_bins = length(bounds) + 1
    #println(num_bins)
    expected_frequencies = fill(total_values / num_bins, num_bins)
    #println(expected_frequencies)
    expected_frequencies = expected_frequencies / total_values
    #println(expected_frequencies)

    return pvalue(ChisqTest(observed_frequencies, expected_frequencies))
end

@everywhere function enroll_optimized(linearcombinations, bounds) # 340ms
    """
    Enrollment function optimized for various input parameters boundsparam. 

    Things that have been done to increase the performance of the original enrollment function 
    - Instead of generating the LinearCombinationSet objects inside of the enrollment function we will instead calculate all these possible results beforehand
    - Currently, there are 8 different combinations to be tested here (n=3 case). Becasue the MSB doesnt really matter for this analysis, we can also lower here the number of input vectors that will be tested. 
    - DO NOT USE pmap, slows down the function by 700ms..
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

println("Starting the brute force computation now.")

# Perform parallel computation with every combination to brute force the result.

#Threads.@threads map(bounds -> begin
#        # Quantize the resulting values and perform the Chi Squared Test
#        pvalue = chi_squared_test(enroll_optimized(linearcombinations, bounds), bounds)
#        # Write the result of p value and the corresponding bounds to the CSV
#        save_to_csv(vcat(pvalue, bounds))
#end, possible_bounds)

Threads.@threads for bounds in possible_bounds
    # Quantize the resulting values and perform the Chi Squared Test
    pvalue = chi_squared_test(enroll_optimized(linearcombinations, bounds), bounds)
    # Write the result of p value and the corresponding bounds to the CSV
    save_to_csv(vcat(pvalue, bounds))
end

# Create a channel to hold the tasks
#task_channel = Channel{Vector{Float64}}(length(possible_bounds))

# Enqueue tasks
#for bounds in possible_bounds
#    put!(task_channel, bounds)
#end

# Function to process tasks from the channel
#function process_tasks(task_channel)
#    for bounds in task_channel
#            bounds = take!(task_channel)
#            pvalue = chi_squared_test(enroll_optimized(linearcombinations, bounds), bounds)
#            save_to_csv(vcat(pvalue, bounds))
#    end
#end

# Launch tasks
#@sync begin
#    for _ in 1:Threads.nthreads()
#        @async process_tasks(task_channel)
#    end
#end

println("Das ist jetzt vorbei hier.")
