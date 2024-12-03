using Distributions

function key_error_rate(n, t, ser, frames=1)
    """Compute the propability that at least one frame of `frames` cannot
    be recovered if `t` symbol/bit errors (occurring with rate `ser`) in
    the each frame of length `n` symbols/bits can be corrected.
    """
    binom_dist = Binomial(n, ser)
    return 1 - cdf(binom_dist, t) ^ frames
end