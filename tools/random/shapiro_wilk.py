import sys
import argparse
import scipy.stats
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("file", type = str, help = "Path of the input file.")

def main():
    args = parser.parse_args();

    with open(args.file, "r") as file:
        values = [float(line) for line in file]

    print(scipy.stats.shapiro(values))

    mean = np.mean(values)
    stdev = np.std(values)

    print(f"mean: {mean}\nstddev: {stdev}")

if __name__ == "__main__":
    main()
