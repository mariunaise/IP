using Pkg
Pkg.activate("env")

using Random
using Distributions

include("julia_code/bach.jl")

function main()
    mean = 0
    std_dev = 1

    dist = Normal(mean, std_dev)
    data = rand(dist, 1000000)

    bach(data)
end

main()