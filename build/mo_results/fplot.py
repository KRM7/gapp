# Script for plotting the fitness values of the population of the NSGA-II and NSGA-III algorithms.

import matplotlib.pyplot as plt
import os

# List of problems to plot the solutions for
problems2 = ["Kursawe", "ZDT1", "ZDT2", "ZDT3", "ZDT4", "ZDT5", "ZDT6"]
problems3 = ["DTLZ1", "DTLZ2", "DTLZ3", "DTLZ4", "DTLZ5", "DTLZ6", "DTLZ7"]


# Read and return data points from file
def readData(name):
    with open(name) as f:
        lines = f.readlines()
        x = [float(line.split()[0]) for line in lines]
        y = [float(line.split()[1]) for line in lines]

    return x, y


# Read and return data points from file
def readData3D(name):
    with open(name) as f:
        lines = f.readlines()
        x = [float(line.split()[0]) for line in lines]
        y = [float(line.split()[1]) for line in lines]
        z = [float(line.split()[2]) for line in lines]

    return x, y, z


# Plot the fitness values of the NSGA-II and III for the same problem
def plotData(name1, name2, title):
    if not os.path.isfile(name1) or not os.path.isfile(name2):
        return

    x1, y1 = readData(name1)
    x2, y2 = readData(name2)

    fig, ax = plt.subplots(1,3)
    ax[0].scatter(x1, y1, c = "red", marker = ".", label = "NSGA-II", linewidths = 1.0)
    ax[0].scatter(x2, y2, c = "blue", marker = ".", label = "NSGA-III", linewidths = 1.0)
    ax[1].scatter(x1, y1, c = "red", marker  = ".", label = "NSGA-II", linewidths = 1.0)
    ax[2].scatter(x2, y2, c = "blue", marker = ".", label = "NSGA-III", linewidths = 1.0)

    for axis in ax:
        axis.set(xlabel = "f1", ylabel = "f2")
        axis.set_title(title)
        axis.invert_xaxis()
        axis.invert_yaxis()
        axis.set_aspect("equal")

    fig.legend(loc = "upper right")
    fig.set_tight_layout(True)
    plt.show()


# Plot the fitness values of the last population
def plotData3D(name1, name2, title):
    if not os.path.isfile(name1) or not os.path.isfile(name2):
        return

    x1, y1, z1 = readData3D(name1)
    x2, y2, z2 = readData3D(name2)

    fig = plt.figure()
    ax = fig.add_subplot(1, 2, 1, projection = "3d")
    ax.scatter(x1, y1, z1, c = "r", marker = "o", label = "NSGA-II", linewidths = 1.5)
    ax.set(xlabel = "f1", ylabel = "f2", zlabel = "f3")
    ax.set_title(title)
    ax.invert_xaxis()
    ax.invert_yaxis()
    ax.invert_zaxis()

    ax = fig.add_subplot(1, 2, 2, projection = "3d")
    ax.scatter(x2, y2, z2, c = "b", marker = "o", label = "NSGA-III", linewidths = 1.5)
    ax.set(xlabel = "f1", ylabel = "f2", zlabel = "f3")
    ax.set_title(title)
    ax.invert_xaxis()
    ax.invert_yaxis()
    ax.invert_zaxis()

    fig.legend(loc = "upper right")
    fig.set_tight_layout(True)
    plt.show()


def main():
    for problem in problems2:
        f1 = "NSGA2_" + problem + "_sols.txt"
        f2 = "NSGA3_" + problem + "_sols.txt"
        plotData(f1, f2, problem + " fitness values")

    for problem in problems3:
        f1 = "NSGA2_" + problem + "_sols.txt"
        f2 = "NSGA3_" + problem + "_sols.txt"
        plotData3D(f1, f2, problem)


if __name__ == "__main__":
    main()