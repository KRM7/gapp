import argparse
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("file", type = str, help = "Path of the input file.")

def main():
    args = parser.parse_args();

    with open(args.file, "r") as file:
        values = [int(line) for line in file]

    plt.hist(values, bins = 200, density = True, label = "Histogram")
    plt.xlabel("x")
    plt.ylabel("p(x)")
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()
