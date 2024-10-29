from dataclasses import dataclass
import numpy as np
import xarray as xr

@dataclass
class LinearCombination:
    weights: list
    value: float 
    

@dataclass
class LinearCombinationSet:
    combination: list

# Function to create a list of weights in +-1
def generate_weights(n):
    weights = np.arange(2**n)
    weights = np.array([[1 if bit == '1' else -1 for bit in np.binary_repr(i, width=n)] for i in weights])
    #weights[:, 0] = 2
    return weights

def generate_cursed_weights(n):
    weights = np.arange(2**n)
    weights = np.array([[2 if bit == '1' else -2 for bit in np.binary_repr(i, width=n)] for i in weights])
    weights[:, 0] = 1
    return weights

# Apply linearcombination takes the input values, weights and number of addends, returns  
def apply_linearcombination(input, weights, n):
    # DataArrays for better multiplication
    input_sliced = xr.DataArray(list(zip(*[iter(input)]*n)), dims=["entry", "values"])
    weights_array = xr.DataArray(weights, dims=["weights", "values"])

    # Expand weights to match the 'entry' dimension of enrolled_sliced
    weights_expanded = weights_array.expand_dims(dim={"entry": input_sliced.sizes["entry"]})

    # Perform element-wise multiplication and sum along the 'values' dimension
    lincombs_raw = (input_sliced * weights_expanded).sum(dim='values')

    # Flatten values down to a single dimension
    lincombs_values = lincombs_raw.values.flatten()
    # Expand weight indicdes 
    weights_indices = np.tile(np.arange(weights_array.sizes['weights']), len(lincombs_raw))

    lincombs = [
        LinearCombination(weights=weights_array.sel(weights=weights_indices[i]).values.tolist(), value=float(lincombs_values[i]))
        for i in range(len(lincombs_values))
    ]

    lincomb_sets = [
        LinearCombinationSet(combination=lincombs[i * weights_array.sizes['weights']:(i+1) * weights_array.sizes['weights']])
        for i in range(len(lincombs_raw))
    ]
    return lincomb_sets, lincombs_values