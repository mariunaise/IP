import numpy as np
import bach, prepare, quantizer_bounds
import imageio
mu, sigma = 0, 1
num_measurements = 1000000
# Enrolled values 
enrolled = np.random.normal(mu, sigma, num_measurements)
weights = prepare.generate_weights(3)
lincombs = prepare.apply_linearcombination(enrolled, weights, 3)
lincombs_values = lincombs[1]
for i in range(4, 100):
    bounds = quantizer_bounds.lazy_with_param(2, lincombs_values, i)
    optimal_combinations = bach.enroll(lincombs[0], bounds)
    # Extract values of linear combinations here again for the plot
    optimal_combinations_raw = [
        combination[0].value
        for combination in optimal_combinations
    ]
    import matplotlib.pyplot as plt

    filenames = []

    # Create histogram plot
    plt.hist(optimal_combinations_raw, bins=1000, alpha=0.75)
    plt.title(f'Histogram for i={i}')
    plt.xlabel('Value')
    plt.ylabel('Frequency')

    # Save the plot as a PNG file
    filename = f'histogram_{i}.png'
    plt.savefig(filename)
    filenames.append(filename)
    plt.close()

    # Create a GIF from the saved PNG files
    if i == 99:
        with imageio.get_writer('histograms.gif', mode='I', duration=0.5) as writer:
            for filename in filenames:
                image = imageio.imread(filename)
                writer.append_data(image)