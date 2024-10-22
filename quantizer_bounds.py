# Set of functions to define the quantizer bounds based on the input values

import numpy as np
import scipy.stats as sct
import scipy as sc

# Default definition based on tile domain equidistant bounds transformed  
def default(m, lincombs_values):
    # Quantizer bounds in the Tilde Domain
    bounds = list(np.linspace(0, 1, 2**m+1))

    print("Bounds in Tilde domain: " , bounds)

    # Bounds real can now be used as constatints for the linear combination to opitmize the input values
    bounds_real= list(map(lambda x: sct.norm.ppf(x, loc=0, scale=lincombs_values.std()), bounds))

    # Remove plus and minus infinity from the list bounds_real
    bounds_real = [bound for bound in bounds_real if not np.isinf(bound)]

    print("Bounds in transformed domain (without infinity): ", bounds_real)

    return bounds_real

# Instead of using the equidistant areas under the PDF of a normal distribution to define the quantizing and optimization bounds for BACH,
# we can try a different definition of these kind of bounds, based on mean and standard deviation of the lincombs_values distribution

def lazy(m, lincombs_values):
    # Define the bounds as a linspace between mu - 3sigma .. mu + 3sigma, lazy halt
    mu = np.average(lincombs_values)
    sigma = np.sqrt(np.std(lincombs_values))

    print(mu, sigma)

    bounds = list(np.linspace(mu - 3*sigma, mu + 3*sigma, 2**m+1))
    print(bounds)
    return bounds
