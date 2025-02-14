### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 21de61ea-ccff-11ef-0c75-496eee3d50da
begin
	using Distributions
	using Random
	using Gadfly
	using LinearAlgebra
	using Statistics
	using DataFrames
	using StatsBase
end

# ╔═╡ 0750906e-480d-4445-91da-5bb3619b838b
# We can test how good the resulting distribution approximates a uniform distribution using the Chi Square Test
using HypothesisTests 

# ╔═╡ 498d4fd1-fa0f-4da4-974a-21092f2860b2
module bach
	include("julia_code/bach_global.jl")
end

# ╔═╡ 6b6f1863-4c23-472d-810a-078862316cef
# Define enrollment values
begin
	dist_mean = 0
	std_dev = 1

	n = 3
	m = 2

	value_len = 1000000
	
	dist = Normal(dist_mean, std_dev)
	data = rand(dist, value_len)
end

# ╔═╡ de1bd394-fa0c-4547-8169-197b372d7ee2
# Define weights 
weights = bach.generate_n_bit_numbers_alpha(n, 1)

# ╔═╡ b30dce14-e736-4ae9-8bee-7f7ae56ff545
linearcombinations = bach.create_linearcombinations(data, weights, n)

# ╔═╡ c53aac8b-6abd-42a1-8e7f-eedf4b1f1dd1
begin 
	iterations = 1
	alpha = 1
end

# ╔═╡ e9bbbfc3-b207-451e-ad28-92261933cba5
# Extract values of linear combinations and do Rayleigh Estimation for the standard deviation 
values = collect(Iterators.flatten(map(set -> begin
	map(combination -> combination.value,set.combination)
end,linearcombinations)))

# ╔═╡ 00a312e7-840b-4f6f-a1af-958f6a531d20
sigma = sqrt(1/2*length(values)*sum(x -> x^2, values))

# ╔═╡ bf8da4fe-6312-4f1d-9b5d-f41b39e83bdb
md"""
-> sigma is not related to the Rayleigh Distribution
"""

# ╔═╡ c0946799-c485-41ef-bb5b-3b97910903db
enrolled = bach.enroll(data, n, m, iterations, alpha, 0.84 * n)

# ╔═╡ f22e78a6-0e69-4052-af7f-4fb9fa7a3cf1
@time bach.enroll(data, n, m, iterations, alpha, 0.84 * n)

# ╔═╡ 84456144-2d49-4e76-ab14-31f970837eaf
plot(x=collect(map(comb -> comb.value, enrolled[1])), Geom.histogram(bincount=1000), Guide.title("Sensible Apprpach"))

# ╔═╡ 1abc1bb5-eb93-406d-8781-161fe7216a37
# Quantize values and plot their histogram:
begin 
	combination_values_1 = map(comb -> comb.value, enrolled[1])
	quantized_indices_1 = searchsortedlast.(Ref(enrolled[2]), combination_values_1)
	plot(x=quantized_indices_1, Geom.histogram(bincount=(length(enrolled[2]) + 1)), Guide.title("Histogram of the quantized values"))
end

# ╔═╡ 8d37dd34-2046-4489-b17e-3e3b167a5a7e
md"""
2.52 seems like the optimal guess for the resulting standard deviation of the result distribution. However, there seems to be no analytical connection to the input distribution.\
This approach of fixing the resulting distribution also does not work for any higher order bit cases higher than 2, since the standrad devation does not have a big enough influence to correct the distribution.
"""

# ╔═╡ 0db60853-92f9-4d87-a718-44beaea71756
dictt = collect(DataFrames.values(countmap(quantized_indices_1)))

# ╔═╡ 8cac6a0a-becf-4e1d-9dde-ff0e1c8e684d
begin
	observed_frequencies = collect(DataFrames.values(countmap(quantized_indices_1)))
	total_values = length(combination_values_1)
	num_bins = length(enrolled[2]) + 1
	expected_frequencies = fill(total_values / num_bins, num_bins)
	expected_frequencies = expected_frequencies / total_values
	#println(expected_frequencies_3)
	#expected_frequencies_3 = [1/4, 1/4, 1/4, 1/4]

	pvalue(ChisqTest(observed_frequencies, expected_frequencies))
end

# ╔═╡ 051e74c7-f7e6-4a10-a678-dd980269819d
md"""
### Helper Data Distribution
"""

# ╔═╡ 4f654b6b-196d-4d32-905e-e227925a8f11
# Function to split the Vector of LinearCombination objects into smaller chunks to which they would be quantized for a better comparison of the distrubtion of helper data
function quantize_objects(objects, bounds)
    bins = [Vector{Any}() for _ in 1:(length(bounds) + 1)]
    for obj in objects
        bin_index = searchsortedlast(bounds, obj.value)
        if obj.value <= bounds[1]
            bin_index = 1  # Place values less than or equal to the smallest bound in the first bin
        elseif obj.value > bounds[end]
            bin_index = length(bins)  # Place values greater than the largest bound in the last bin
        else
            bin_index += 1  # Adjust bin index for values within the bounds
        end
        push!(bins[bin_index], obj)
    end
    return bins
end

# ╔═╡ 16387326-9947-4246-866f-5be13b7e158f
quantized_objects = quantize_objects(enrolled[1], enrolled[2])

# ╔═╡ d6427f24-5bc9-4c09-afa7-076e18733dc4
overall_helperdata = map(v -> v.helperdatabits ,enrolled[1])

# ╔═╡ c4d2825c-b925-4fb0-839a-3c9a8bc1b1bc
hdocs = map(bin -> begin 
	map(v -> v.helperdatabits ,bin)
end,quantized_objects)

# ╔═╡ f280afa1-fc05-4dae-a741-7bfa3e1819b3
occurs = map(ocs -> countmap(ocs), hdocs)

# ╔═╡ d5b76d16-f26a-47f8-b5cd-78045c2a5350
counts = countmap(overall_helperdata)

# ╔═╡ d1e5d390-f7b7-4a61-baf0-159a1e58d9e6
begin
	# Convert the dictionary to a DataFrame for plotting
	bool_vectors = [string(k) for k in DataFrames.keys(counts)]
	counts_values = collect(DataFrames.values(counts))
	df = DataFrame(BoolVector = bool_vectors, Count = counts_values)

	# Create a bar plot
	p = plot(df, x=:BoolVector, y=:Count, Geom.bar, 
         Guide.xlabel("Bool Vector"), Guide.ylabel("Count"), 
         Guide.title("Occurrences of Bool Vectors"))
end

# ╔═╡ 31df2498-80dc-412c-8065-838038ac76e6
begin
combined_data = DataFrame(BoolVector = String[], Count = Int[], Source = String[])

for (i, dict) in enumerate(occurs)
    bool_vectors = [string(k) for k in DataFrames.keys(dict)]
    counts_values = collect(DataFrames.values(dict))
    source_label = fill("Symbol $i", length(bool_vectors))
    combined_data = vcat(combined_data, DataFrame(BoolVector = bool_vectors, Count = counts_values, Source = source_label))
end

# Create a bar plot with different colors for each source
	plot(combined_data, x=:BoolVector, y=:Count, color=:Source, Geom.bar, 
         Guide.xlabel("Bool Vector"), Guide.ylabel("Count"), 
         Guide.title("Helper Data occurrences by quantized bits"))
end

# ╔═╡ d7b973e7-cbcd-411e-9c3d-2cab4fa63326
md"""
# YES
"""

# ╔═╡ 5e43e2ac-b040-4ea0-9316-d8f18ee2aa9e
md"""
## Investitating the distribution parameters of the resulting distribution
"""

# ╔═╡ 5185b144-5294-4e83-89dc-2fa864bc08d1
# Standard deviation of the resulting values after filtering
sigma_new = Statistics.std(map(comb -> comb.value, enrolled[1]))

# ╔═╡ ef73c678-db04-4672-86dc-de6f34a8e382
enrolled_corrected = bach.enroll(data, n, m, iterations, alpha, sigma_new)

# ╔═╡ c5eed418-00fb-47ad-9f90-ba82dd523d15
# Quantize values and plot their histogram:
begin 
	combination_values_2 = map(comb -> comb.value, enrolled_corrected[1])
	quantized_indices_2 = searchsortedlast.(Ref(enrolled_corrected[2]), combination_values_2)
	plot(x=quantized_indices_2, Geom.histogram(bincount=(length(enrolled_corrected[2]) + 1)), Guide.title("Histogram of the quantized values"))
end

# ╔═╡ 8032e836-b549-461f-a2f8-aa6087952aca
begin
	observed_frequencies_2 = collect(DataFrames.values(countmap(quantized_indices_2)))
	total_values_2 = length(combination_values_2)
	num_bins_2 = length(enrolled[2]) + 1
	expected_frequencies_2 = fill(total_values / num_bins, num_bins)

	chi_square_stat_2 = sum((observed_frequencies_2 .- expected_frequencies_2) .^2 ./ expected_frequencies_2)

	p_value_2 = 1 - cdf(Chisq(num_bins_2 -1), chi_square_stat_2)
end

# ╔═╡ 0c001087-fb4f-41ca-8c2b-2de50a168078
# Standard deviation of the resulting values after filtering again
sigma_new2 = Statistics.std(map(comb -> comb.value, enrolled_corrected[1]))

# ╔═╡ 2bae9b27-344b-49f0-bc70-f2b038aba29e
enrolled_corrected2 = bach.enroll(data, n, m, iterations, alpha, sigma_new)

# ╔═╡ b65b20c5-1c11-4aba-8460-4ec5496c541f
# Quantize values and plot their histogram:
begin 
	combination_values_3 = map(comb -> comb.value, enrolled_corrected2[1])
	quantized_indices_3 = searchsortedlast.(Ref(enrolled_corrected2[2]), combination_values_3)
	plot(x=quantized_indices_3, Geom.histogram(bincount=(length(enrolled_corrected2[2]) + 1)), Guide.title("Histogram of the quantized values"))
end

# ╔═╡ 84015914-9f1c-4bb0-9226-bb1658f67c21


# ╔═╡ 69597adb-c019-4158-8363-47824cc0ca9f
begin
	observed_frequencies_3 = collect(DataFrames.values(countmap(quantized_indices_3)))
	total_values_3 = length(combination_values_3)
	num_bins_3 = length(enrolled[2]) + 1
	expected_frequencies_3 = fill(total_values_3 / num_bins_3, num_bins_3)
	expected_frequencies_3 = expected_frequencies_3 / total_values
	#println(expected_frequencies_3)
	#expected_frequencies_3 = [1/4, 1/4, 1/4, 1/4]

	ChisqTest(observed_frequencies_3, expected_frequencies_3)
end

# ╔═╡ 4c778ed5-1dff-4d89-a7e7-b880d43fb3c2
# Standard deviation of the resulting values after filtering again
Statistics.std(map(comb -> comb.value, enrolled_corrected[1]))

# ╔═╡ 5364b3d0-2695-43dd-b6bc-5fa228bd7a70
md"""
-> Iterative approach to find the sigma does not find the correct solution
"""

# ╔═╡ 223b2fed-22bf-4723-89d7-e8671dd419c3
md"""
# Finding bounds based on the recursive approach and using them in the global variant
"""

# ╔═╡ ea615cf1-ae76-4b06-b89b-be5b0d7319fb
module bach_recursive 
	include("julia_code/bach.jl")
end

# ╔═╡ 27c75d71-6fed-4d23-b924-263f8a1f3a96
function number_sequence(i)
	1 / 2^i
end

# ╔═╡ d8fb0b7f-77e7-4cce-94e0-bdd2be432fc1
bounds = bach_recursive.bach(data, n, m, number_sequence)[2]

# ╔═╡ e21faef2-675a-4b6d-9447-2201a6886f70
bounds_rounded = round.(bounds, digits=2)

# ╔═╡ adfdf1fb-7454-4f3d-8907-7e5a8b47e730
enrolled_with_recursive_bounds = bach.enroll(data, n, m, iterations, alpha, 3.2, bounds_rounded)[1]

# ╔═╡ 65b85096-e0a6-4057-906f-59afd9354d18
plot(x=collect(map(comb -> comb.value, enrolled_with_recursive_bounds)), Geom.histogram(bincount=1000), Guide.title("Sensible Apprpach"))

# ╔═╡ 060c016f-c9c0-4b0f-832d-093199ba61b6
md"""
Take a look at the quantized distribution of bits
"""

# ╔═╡ c97ef68d-0286-4691-9a2c-97dc9a63bf20
# Quantize values and plot their histogram:
begin 
	combination_values = map(comb -> comb.value, enrolled_with_recursive_bounds)
	quantized_indices = searchsortedlast.(Ref(bounds_rounded), combination_values)
	plot(x=quantized_indices, Geom.histogram(bincount=(length(bounds_rounded) + 1)), Guide.title("Histogram of the quantized values"))
end

# ╔═╡ 37b70444-ca73-424b-81dd-4a2bfb6d2e05
# Standard deviation of the resulting values after filtering again
Statistics.std(map(comb -> comb.value, enrolled_with_recursive_bounds))

# ╔═╡ cb9a75e9-76d3-4239-9553-72a5eeeeed95
md"""
As we can see, this approach does not instantly yield optmial bounds
"""

# ╔═╡ 39da5c79-8665-4090-91b3-72f0cea997e2
md"""
# BER performance analysis
"""

# ╔═╡ 12d8a404-8d90-442f-8c74-e3ca7853070d
# Define distribution parameters for reconstruction values 
begin 
	error_dist = Normal(0, 0.01)
	error_values = rand(error_dist, value_len)

	recon_data = data .+ error_values
end

# ╔═╡ bd7fea5f-95a2-4ab9-ba0e-f2abd78af3db
found_weights = map(comb -> comb.weights, enrolled[1])

# ╔═╡ 034cab50-5a94-4b9c-8fb1-7b467b91d8fa
found_bounds = enrolled[2]

# ╔═╡ 616638c1-b4f6-4a32-af3f-4ff5ffa0b925
enrolled_values = map(comb -> comb.value, enrolled[1])

# ╔═╡ ec6daa43-58e6-4de2-a77f-b3d126c20dcb
reconstructed = bach.reconstruct(recon_data,n, found_weights)

# ╔═╡ d6832adf-847a-4822-abf2-7d4c51171e74
md"""
### Compare enrolled and reconstructed quantized codewords
"""

# ╔═╡ 4f00c0ca-9c92-4139-bb65-633eebf996c9
enrolled_codeword = searchsortedlast.(Ref(found_bounds), enrolled_values)

# ╔═╡ a83a36a7-1517-4700-b5f9-1ef853dfcc5e
reconstructed_codeword = searchsortedlast.(Ref(found_bounds), reconstructed)

# ╔═╡ 33935c58-c81e-4b5a-bcec-3aafa127cfe9
md"""
### Calculate the SER of this BACH optimization
"""

# ╔═╡ ef6407e5-3980-4b8f-bff0-1333479cfe91
error_count = count(part -> begin 
	part[1] != part[2]
end,zip(enrolled_codeword, reconstructed_codeword))

# ╔═╡ 2ef511e7-7ecd-4c80-b782-ef4d82429737
error_rate = error_count / length(enrolled_codeword)

# ╔═╡ e73762de-ac03-4881-8c0a-7505ab9ce01e
md"""
### Comparison, if we would do nothing 
"""

# ╔═╡ 12af5893-3506-4e12-917f-c39da6843ce8
nothing_codeword = searchsortedlast.(Ref(found_bounds), data)

# ╔═╡ a1200828-a0ae-4c41-ada3-03129b94b6b1
nothing_reconstructed = searchsortedlast.(Ref(found_bounds), recon_data)

# ╔═╡ ede6b117-8107-425e-8aeb-d9cc9d862f93
nothing_error_count = count(part -> begin 
	part[1] != part[2]
end,zip(nothing_codeword, nothing_reconstructed))

# ╔═╡ 8659f8ff-029b-406f-b7d4-a09191ffbb3b
nothing_error_rate = nothing_error_count / length(data)

# ╔═╡ 3243be6f-a646-4cc0-a803-d4a5b0f1060b
md"""
# Thank god its better
"""

# ╔═╡ b6c71800-acd3-49d5-817a-b69f42b3e32b
md"""
Visual of the found bound solution
"""

# ╔═╡ d6b2ed89-6415-4207-8552-0c2c2ef08c7a
begin 
	a = 2.71
	b = 1.67
	c = 0.81
	solutions = [-a, -b, -c, 0.0, c, b, a]
end

# ╔═╡ 3ff58c0c-31c7-4990-bffc-9269e0b2dfb2
enrolled_higher_order = bach.enroll(data, n, 3, iterations, alpha, 3.2, solutions)[1]

# ╔═╡ 184fa6e0-caa3-46bc-96b6-0ccaca7bd040
plot(x=collect(map(comb -> comb.value, enrolled_higher_order)), Geom.histogram(bincount=1000), Guide.title("Sensible Apprpach"))

# ╔═╡ ee3ee264-1a63-4e17-9830-916335154f65
# Quantize values and plot their histogram:
begin 
	combination_values_5 = map(comb -> comb.value, enrolled_higher_order)
	quantized_indices_5 = searchsortedlast.(Ref(solutions), combination_values_5)
	plot(x=quantized_indices_5, Geom.histogram(bincount=(length(solutions) + 1)), Guide.title("Histogram of the quantized values"))
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
Gadfly = "c91e804a-d5a3-530f-b6f0-dfbca275c004"
HypothesisTests = "09f84164-cd44-5f33-b23f-e6b0d136a0d5"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
DataFrames = "~1.7.0"
Distributions = "~0.25.114"
Gadfly = "~1.4.0"
HypothesisTests = "~0.11.3"
Statistics = "~1.11.1"
StatsBase = "~0.33.21"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.3"
manifest_format = "2.0"
project_hash = "23eeef4b8307fe126926cc75dd1986bef3a98047"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "0ba8f4c1f06707985ffb4804fdad1bf97b233897"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.41"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "50c3c56a52972d78e8be9fd135bfb91c9574c140"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.1.1"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.CategoricalArrays]]
deps = ["DataAPI", "Future", "Missings", "Printf", "Requires", "Statistics", "Unicode"]
git-tree-sha1 = "1568b28f91293458345dabba6a5ea3f183250a61"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.10.8"

    [deps.CategoricalArrays.extensions]
    CategoricalArraysJSONExt = "JSON"
    CategoricalArraysRecipesBaseExt = "RecipesBase"
    CategoricalArraysSentinelArraysExt = "SentinelArrays"
    CategoricalArraysStructTypesExt = "StructTypes"

    [deps.CategoricalArrays.weakdeps]
    JSON = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SentinelArrays = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
    StructTypes = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "3e4b134270b372f2ed4d4d0e936aabaefc1802bc"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Compose]]
deps = ["Base64", "Colors", "DataStructures", "Dates", "IterTools", "JSON", "LinearAlgebra", "Measures", "Printf", "Random", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "bf6570a34c850f99407b494757f5d7ad233a7257"
uuid = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
version = "0.9.5"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConstructionBase]]
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CoupledFields]]
deps = ["LinearAlgebra", "Statistics", "StatsBase"]
git-tree-sha1 = "6c9671364c68c1158ac2524ac881536195b7e7bc"
uuid = "7ad07ef1-bdf2-5661-9d2b-286fd4296dac"
version = "0.2.0"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "fb61b4812c49343d7ef0b533ba982c46021938a6"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.7.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.12"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "9d9e93d19c912ee6f0f3543af0d8839079dbd0d7"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.114"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "4820348781ae578893311153d69049a93d05f39d"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4d81ed14783ec49ce9f2e168208a12ce1815aa25"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.Gadfly]]
deps = ["Base64", "CategoricalArrays", "Colors", "Compose", "Contour", "CoupledFields", "DataAPI", "DataStructures", "Dates", "Distributions", "DocStringExtensions", "Hexagons", "IndirectArrays", "IterTools", "JSON", "Juno", "KernelDensity", "LinearAlgebra", "Loess", "Measures", "Printf", "REPL", "Random", "Requires", "Showoff", "Statistics"]
git-tree-sha1 = "d546e18920e28505e9856e1dfc36cff066907c71"
uuid = "c91e804a-d5a3-530f-b6f0-dfbca275c004"
version = "1.4.0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.Hexagons]]
deps = ["Test"]
git-tree-sha1 = "de4a6f9e7c4710ced6838ca906f81905f7385fd6"
uuid = "a1b4810d-1bce-5fbd-ac56-80944d57a21f"
version = "0.2.0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "b1c2585431c382e3fe5805874bda6aea90a95de9"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.25"

[[deps.HypothesisTests]]
deps = ["Combinatorics", "Distributions", "LinearAlgebra", "Printf", "Random", "Rmath", "Roots", "Statistics", "StatsAPI", "StatsBase"]
git-tree-sha1 = "6c3ce99fdbaf680aa6716f4b919c19e902d67c9c"
uuid = "09f84164-cd44-5f33-b23f-e6b0d136a0d5"
version = "0.11.3"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.InlineStrings]]
git-tree-sha1 = "45521d31238e87ee9f9732561bfee12d4eebd52d"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.2"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "10bd689145d2c3b2a9844005d01087cc1194e79e"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.2.1+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"

    [deps.Interpolations.extensions]
    InterpolationsUnitfulExt = "Unitful"

    [deps.Interpolations.weakdeps]
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "be3dc50a92e5a386872a493a10050136d4703f9b"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.6.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.Juno]]
deps = ["Base64", "Logging", "Media", "Profile"]
git-tree-sha1 = "07cb43290a840908a771552911a6274bc6c072c7"
uuid = "e5e0dc1b-0480-54bc-9374-aad01c23163d"
version = "0.8.4"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "7d703202e65efa1369de1279c162b915e245eed1"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.9"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.Loess]]
deps = ["Distances", "LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "f749e7351f120b3566e5923fefdf8e52ba5ec7f9"
uuid = "4345ca2d-374a-55d4-8d30-97f9976e7612"
version = "0.6.4"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "f046ccd0c6db2832a9f639e2c669c6fe867e5f4f"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.2.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Media]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "75a54abd10709c01f1b86b84ec225d26e840ed58"
uuid = "e89f7d12-3494-54d1-8411-f7d8b9ae1f27"
version = "0.5.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OffsetArrays]]
git-tree-sha1 = "39d000d9c33706b8364817d8894fae1548f40295"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.14.2"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "12f1439c4f986bb868acda6ea33ebc78e19b95ad"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.7.0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "1101cd475833706e4d0e7b122218257178f48f34"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Profile]]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "77a42d78b6a92df47ab37e177b2deac405e1c88f"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.2.1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "cda3b045cf9ef07a08ad46731f5a3165e56cf3da"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.1"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "852bd0f55565a9e973fcfee83a84413270224dc4"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.8.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "f233e0a3de30a6eed170b8e1be0440f732fdf456"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.2.4"

    [deps.Roots.extensions]
    RootsChainRulesCoreExt = "ChainRulesCore"
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Roots.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "712fb0231ee6f9120e005ccd56297abbc053e7e0"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.8"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "64cca0c26b4f31ba18f13f6c12af7c85f478cfde"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.5.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "777657803913ffc7e8cc20f0fd04b634f871af8f"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.8"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "b423576adc27097764a90e163157bcfc9acf0f46"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.2"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a6b1675a536c5ad1a60e5a5153e1fee12eb146e3"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.0"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7d0ea0f4895ef2f5cb83645fa689e52cb55cf493"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2021.12.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╠═21de61ea-ccff-11ef-0c75-496eee3d50da
# ╠═498d4fd1-fa0f-4da4-974a-21092f2860b2
# ╠═6b6f1863-4c23-472d-810a-078862316cef
# ╠═de1bd394-fa0c-4547-8169-197b372d7ee2
# ╠═b30dce14-e736-4ae9-8bee-7f7ae56ff545
# ╠═c53aac8b-6abd-42a1-8e7f-eedf4b1f1dd1
# ╠═e9bbbfc3-b207-451e-ad28-92261933cba5
# ╠═00a312e7-840b-4f6f-a1af-958f6a531d20
# ╟─bf8da4fe-6312-4f1d-9b5d-f41b39e83bdb
# ╠═c0946799-c485-41ef-bb5b-3b97910903db
# ╠═f22e78a6-0e69-4052-af7f-4fb9fa7a3cf1
# ╠═84456144-2d49-4e76-ab14-31f970837eaf
# ╠═1abc1bb5-eb93-406d-8781-161fe7216a37
# ╟─8d37dd34-2046-4489-b17e-3e3b167a5a7e
# ╠═0750906e-480d-4445-91da-5bb3619b838b
# ╠═0db60853-92f9-4d87-a718-44beaea71756
# ╠═8cac6a0a-becf-4e1d-9dde-ff0e1c8e684d
# ╟─051e74c7-f7e6-4a10-a678-dd980269819d
# ╟─4f654b6b-196d-4d32-905e-e227925a8f11
# ╠═16387326-9947-4246-866f-5be13b7e158f
# ╠═d6427f24-5bc9-4c09-afa7-076e18733dc4
# ╠═c4d2825c-b925-4fb0-839a-3c9a8bc1b1bc
# ╠═f280afa1-fc05-4dae-a741-7bfa3e1819b3
# ╠═d5b76d16-f26a-47f8-b5cd-78045c2a5350
# ╠═d1e5d390-f7b7-4a61-baf0-159a1e58d9e6
# ╟─31df2498-80dc-412c-8065-838038ac76e6
# ╟─d7b973e7-cbcd-411e-9c3d-2cab4fa63326
# ╟─5e43e2ac-b040-4ea0-9316-d8f18ee2aa9e
# ╠═5185b144-5294-4e83-89dc-2fa864bc08d1
# ╠═ef73c678-db04-4672-86dc-de6f34a8e382
# ╠═c5eed418-00fb-47ad-9f90-ba82dd523d15
# ╠═8032e836-b549-461f-a2f8-aa6087952aca
# ╠═0c001087-fb4f-41ca-8c2b-2de50a168078
# ╠═2bae9b27-344b-49f0-bc70-f2b038aba29e
# ╠═b65b20c5-1c11-4aba-8460-4ec5496c541f
# ╠═84015914-9f1c-4bb0-9226-bb1658f67c21
# ╠═69597adb-c019-4158-8363-47824cc0ca9f
# ╠═4c778ed5-1dff-4d89-a7e7-b880d43fb3c2
# ╟─5364b3d0-2695-43dd-b6bc-5fa228bd7a70
# ╟─223b2fed-22bf-4723-89d7-e8671dd419c3
# ╠═ea615cf1-ae76-4b06-b89b-be5b0d7319fb
# ╠═27c75d71-6fed-4d23-b924-263f8a1f3a96
# ╠═d8fb0b7f-77e7-4cce-94e0-bdd2be432fc1
# ╠═e21faef2-675a-4b6d-9447-2201a6886f70
# ╠═adfdf1fb-7454-4f3d-8907-7e5a8b47e730
# ╠═65b85096-e0a6-4057-906f-59afd9354d18
# ╠═060c016f-c9c0-4b0f-832d-093199ba61b6
# ╠═c97ef68d-0286-4691-9a2c-97dc9a63bf20
# ╟─37b70444-ca73-424b-81dd-4a2bfb6d2e05
# ╟─cb9a75e9-76d3-4239-9553-72a5eeeeed95
# ╟─39da5c79-8665-4090-91b3-72f0cea997e2
# ╠═12d8a404-8d90-442f-8c74-e3ca7853070d
# ╠═bd7fea5f-95a2-4ab9-ba0e-f2abd78af3db
# ╠═034cab50-5a94-4b9c-8fb1-7b467b91d8fa
# ╠═616638c1-b4f6-4a32-af3f-4ff5ffa0b925
# ╠═ec6daa43-58e6-4de2-a77f-b3d126c20dcb
# ╟─d6832adf-847a-4822-abf2-7d4c51171e74
# ╠═4f00c0ca-9c92-4139-bb65-633eebf996c9
# ╠═a83a36a7-1517-4700-b5f9-1ef853dfcc5e
# ╟─33935c58-c81e-4b5a-bcec-3aafa127cfe9
# ╠═ef6407e5-3980-4b8f-bff0-1333479cfe91
# ╠═2ef511e7-7ecd-4c80-b782-ef4d82429737
# ╟─e73762de-ac03-4881-8c0a-7505ab9ce01e
# ╠═12af5893-3506-4e12-917f-c39da6843ce8
# ╠═a1200828-a0ae-4c41-ada3-03129b94b6b1
# ╠═ede6b117-8107-425e-8aeb-d9cc9d862f93
# ╠═8659f8ff-029b-406f-b7d4-a09191ffbb3b
# ╟─3243be6f-a646-4cc0-a803-d4a5b0f1060b
# ╠═b6c71800-acd3-49d5-817a-b69f42b3e32b
# ╠═d6b2ed89-6415-4207-8552-0c2c2ef08c7a
# ╠═3ff58c0c-31c7-4990-bffc-9269e0b2dfb2
# ╠═184fa6e0-caa3-46bc-96b6-0ccaca7bd040
# ╠═ee3ee264-1a63-4e17-9830-916335154f65
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
