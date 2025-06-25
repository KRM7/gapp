import argparse
import scipy.stats
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("file", type = str, help = "Path of the input file.")
parser.add_argument("dist", type = str, help = "The expected continuous distribution of the values in the input file.",
                    choices = ["uniform", "normal", "exponential"])

def main():
    args = parser.parse_args();

    with open(args.file, "r") as file:
        values = [float(line) for line in file]

    expected = args.dist

    mean = np.mean(values)
    sdev = np.std(values)

    if expected == "normal":
        expected = "norm"
        values = [(val - mean) / sdev for val in values]

    elif expected == "exponential":
        expected = "expon"
        values = [val / mean for val in values]

    print(scipy.stats.kstest(values, expected))
    print(f"mean: {mean}\nstddev: {sdev}")

if __name__ == "__main__":
    main()
