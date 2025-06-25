import argparse
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("file", type = str, help = "Path of the input file.")

def main():
    args = parser.parse_args();

    with open(args.file, "r") as file:
        values = [float(line) for line in file]

    kde = scipy.stats.gaussian_kde(values)

    x_grid = np.linspace(min(values), max(values), 200)
    pdf_values = kde(x_grid)

    plt.plot(x_grid, pdf_values, label = "Estimated PDF")
    plt.hist(values, bins = 200, density = True, alpha = 0.4, label = "Histogram")
    plt.xlabel("x")
    plt.ylabel("p(x)")
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()
