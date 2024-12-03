using LinearAlgebra
using Random

function poisson_binomial_cdf(x, p)
    n = size(p, 2)
    if x < 0 
        trow(ArgumentError("x must be non-negative."))
    end
    if x > n 
        throw(ArgumentError("x must not be larger than n."))
    end

    C = exp(2im * pi / (n +1))
    nrange = 1:n 

    result = (x + 1) / (n + 1) + 1 / (n + 1) * sum(
        ((1 .- C .^ (-nrange * (x + 1))) ./ (1 .- C .^ (-nrange)))' .* 
        prod(p * C .^ nrange' .+ (1 .- p), dims=2),
        dims=2
    )

    return real(result)
end

function nonuniform_key_error_rate(n, t, sers, frames, samples, memory_limit)
    kers = Array{Float64}(undef, samples)
    
    # Split kers into batches
    num_batches = min(n * frames * samples ÷ memory_limit + 1, samples)
    batch_size = ceil(Int, samples / num_batches)
    batches = [kers[i:min(i + batch_size - 1, samples)] for i in 1:batch_size:samples]

    for batch in batches
        batch_samples = length(batch)
        # Generate random indices
        random_indices = rand(1:size(sers, 1), frames * batch_samples * n)
        reshaped_sers = reshape(sers[random_indices], frames * batch_samples, n)
        
        # Compute the Poisson-Binomial CDF
        cdf_values = poisson_binomial_cdf(t, reshaped_sers)
        reshaped_cdf = reshape(cdf_values, frames, batch_samples)
        
        # Compute the product and update the batch
        batch .= 1 .- prod(reshaped_cdf, dims=1)
    end
    
    return kers
end