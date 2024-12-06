using Random
using Statistics

function repetition_bers(enrollment_data, reconstruction_data, rep_bits, mc_samples, random_references, seed, memory_limit)
    if size(enrollment_data, 1) != size(reconstruction_data, 1)
        throw(ArgumentError("Enrollment and reconstruction data must have the same number of cells."))
    end 
    
    cell_count = size(enrollment_data, 1)

    if cell_count < rep_bits + 1 
        throw(ArgumentError("Not enough cells for chosen repetition code."))
    end

    if rep_bits % 2 == 0
        throw(ArgumentError("Even-sized repetition codes are not supported."))
    end 
    
    if mc_samples * cell_mount > memory_limit
        bers = Array{Float64}(undef, mc_samples)
        chunk_size = Int(div(memory_limit, cell_count))
        
        for i in 1:chunk_size:mc_samples 
            bers[i:min(i + chunk_size -1, mc_samples)] .= repetition_bers(
                enrollment_data, reconstruction_data, rep_bits,
                length(bers[i:min(i + chunk_size - 1, mc_samples)]),
                random_references, 
                seed + i - 1 if seed != nothing else nothing,
                memory_limit
            )
        end
        return bers
    end

    cell_indices = repeat(reshape(collect(0:cell_count-1), 1, :), mc_samples, 1)

    rng = MersenneTwister(seed)

    cell_indices = permutedims(cell_indices, [2, 1])[:, 1:rep_bits + 1]
    
    enrollment_data = enrollment_data[cell_indices, :]
    reconstruction_data = reconstruction_data[cell_indices, :]

    enrollment_means = mean(enrollment_data, dims=3)

    if !random_references 
       
        reference_sort = sortperm(enrollment_means, dims=2)

        cell_permutation = collect(0:rep_bits)
        cell_permutation[1:(rep_bits) / 2 +1] .= 1 
        cell_permutation[1] = (rep_bits + 1) / 2 

        refrence_sort = reference_sort[:, cell_permutation .+ 1]
    
        enrollment_means = enrollment_means[:, reference_sort]
        reconstruction_data = reconstruction_data[:, reference_sort, :]
    end

    # [MC sample, cell]
    enrolment_bits = enrolment_means[:, 2:end] .> enrolment_means[:, 1]

    # [MC sample, cell, sample]
    reconstruction_bits = reconstruction_data[:, 2:end, :] .> reconstruction_data[:, 1, :]

    # [MC sample, cell, sample]
    bit_errors = reconstruction_bits .!= enrolment_bits[:, :, :]

    # Calculate bit errors after repetition decoder and average over measurement samples
    return mean(sum(bit_errors, dims=2) .> rep_bits ÷ 2, dims=2) 
end
