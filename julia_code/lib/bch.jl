include("functions.jl")

using StatsBase

BCH_CODES = [
    (1, 1, 0),
    (7, 4, 1),
    (7, 1, 3),
    (15, 11, 1),
    (15, 7, 2),
    (15, 5, 3),
    (15, 1, 7),
    (31, 26, 1),
    (31, 21, 2),
    (31, 16, 3),
    (31, 11, 5),
    (31, 6, 7),
    (31, 1, 15),
    (63, 57, 1),
    (63, 51, 2),
    (63, 45, 3),
    (63, 39, 4),
    (63, 36, 5),
    (63, 30, 6),
    (63, 24, 7),
    (63, 18, 10),
    (63, 16, 11),
    (63, 10, 13),
    (63, 7, 15),
    (63, 1, 31),
    (127, 120, 1),
    (127, 113, 2),
    (127, 106, 3),
    (127, 99, 4),
    (127, 92, 5),
    (127, 85, 6),
    (127, 78, 7),
    (127, 71, 9),
    (127, 64, 10),
    (127, 57, 11),
    (127, 50, 13),
    (127, 43, 14),
    (127, 36, 15),
    (127, 29, 21),
    (127, 22, 23),
    (127, 15, 27),
    (127, 8, 31),
    (127, 1, 63),
    (255, 247, 1),
    (255, 239, 2),
    (255, 231, 3),
    (255, 223, 4),
    (255, 215, 5),
    (255, 207, 6),
    (255, 199, 7),
    (255, 191, 8),
    (255, 187, 9),
    (255, 179, 10),
    (255, 171, 11),
    (255, 163, 12),
    (255, 155, 13),
    (255, 147, 14),
    (255, 139, 15),
    (255, 131, 18),
    (255, 123, 19),
    (255, 115, 21),
    (255, 107, 22),
    (255, 99, 23),
    (255, 91, 25),
    (255, 87, 26),
    (255, 79, 27),
    (255, 71, 29),
    (255, 63, 30),
    (255, 55, 31),
    (255, 47, 42),
    (255, 45, 43),
    (255, 37, 45),
    (255, 29, 47),
    (255, 21, 55),
    (255, 13, 59),
    (255, 9, 63),
    (255, 1, 127)
]

for (n, k, t) in BCH_CODES
    @assert t <= (n - k + 1) ÷ 2
end

function find_best(bers, rep_bits, key_size, ker_limit, rating, percentile, safety_factor, extra_codes, use_nonuniform, samples)
    """Find suitable BCH codes and return the best ones, i.e. the codes using the fewest memristors to achieve 
    the specifications.
    """
    codes = []
    for (n, k, t) in BCH_CODES
        frames = ceil(Int, key_size / k)
        memristors = n * frames * (rep_bits + 1)
        append!(codes, Dict(
            "code" => (n, k, t),
            "frames" => frames, 
            "memristors" => memristors, 
            "rep_bits" => rep_bits
        ))
    end
    sort!(codes, by = x -> x["memristors"])
    results = []
    for c in codes 
        n, _, t = c["code"]
        if use_nonuniform
            kers = nonuniform_key_error_rate(n, t, bers, frames, samples)
        else
            kers = key_error_rate(n, t, bers, frames)
        end
        mean_ker = mean(kers)
        max_ker = max(kers)
        percentile_ker = percentile(kers, percentile)
        if rating == "mean"
            rated_ker = mean_ker
        elseif rating == "percentile"
            rated_ker = percentile_ker
        else rating == "max"
            rated_ker = max_ker
        end
        if rated_ker * safety_factor < ker_limit 
            c["mean_ker"] = mean_ker
            c["percentile_ker"] = percentile_ker
            c["max_ker"] = max_ker
            push!(results, c)
            extra_codes -= 1
            if extra_codes < 0
                break
            end
        end
    end
    return results
end
