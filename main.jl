using Pkg
Pkg.activate("env")

using Random
using Distributions
ENV["GDK_BACKEND"] = "x11"
using Gadfly

include("julia_code/bach.jl")
include("julia_code/helperfunctions.jl")

function main()

    # Define the parameters of the input probability distribution
    mean = 0
    std_dev = 1

    # Amount of addends in the linear combination
    n = 3

    # Define the number of bits we want to extract
    m = 2

    dist = Normal(mean, std_dev)
    data = rand(dist, 1000000)

    # Call the BACH function to do the magic
    bach( data,n, m)
end

main()