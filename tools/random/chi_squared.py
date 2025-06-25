import sys
import math
import argparse
import collections
import numpy as np
import scipy.stats

parser = argparse.ArgumentParser()
parser.add_argument("file", type = str, help = "Path of the input file.")
parser.add_argument("dist", type = str, help = "The expected discrete distribution of the values in the input file.",
                    choices = ["uniform", "binomial", "poisson"])

def main():
    args = parser.parse_args();

    with open(args.file, "r") as file:
        values = [int(line) for line in file]

    counts = collections.Counter(values)
    counts = {k: v for k, v in counts.items() if v >= 5}

    if args.dist == "uniform":
        print(scipy.stats.chisquare(list(counts.values())))
        print(f"min: {np.min(values)}\nmax: {np.max(values)}")
        sys.exit(0)

    observed = list(counts.values())
    observed = [n / sum(observed) for n in observed]

    mean = np.mean(values)

    if args.dist == "poisson":
        expected = {k: scipy.stats.poisson.pmf(k, mean) for k, v in counts.items() }
        expected = list(expected.values())
        expected = [n / sum(expected) for n in expected]
        print(scipy.stats.chisquare(observed, expected))
        print(f"mean: {mean}")
        sys.exit(0)

    var = np.var(values)
    p = 1 - var / mean
    n = int(round(mean / p))

    if args.dist == "binomial":
        expected = {k: scipy.stats.binom.pmf(k, n, p) for k, v in counts.items() }
        expected = list(expected.values())
        expected = [n / sum(expected) for n in expected]
        print(scipy.stats.chisquare(observed, expected))
        print(f"mean: {mean}\nn: {n}\np: {p}")
        sys.exit(0)

if __name__ == "__main__":
    main()
