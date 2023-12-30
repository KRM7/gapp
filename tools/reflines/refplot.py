# Script for plotting reference lines

import matplotlib.pyplot as plt
import os
import sys


def readData2D(name):
    with open(name) as f:
        lines = f.readlines()
        x = [float(line.split()[0]) for line in lines]
        y = [float(line.split()[1]) for line in lines]

    return x, y


def readData3D(name):
    with open(name) as f:
        lines = f.readlines()
        x = [float(line.split()[0]) for line in lines]
        y = [float(line.split()[1]) for line in lines]
        z = [float(line.split()[2]) for line in lines]

    return x, y, z


def plotData(fname):
    
    fig = plt.figure()

    x2, y2 = readData2D(fname + "_2D.txt")

    ax = fig.add_subplot(1, 2, 1)
    ax.scatter(x2, y2, c = "b", marker = "x", label = "2D", linewidths = 1.5)
    ax.set(xlabel = "x", ylabel = "y")
    ax.set_aspect("equal")

    x3, y3, z3 = readData3D(fname + "_3D.txt")

    ax = fig.add_subplot(1, 2, 2, projection = "3d")
    ax.scatter(x3, y3, z3, c = "b", marker = "x", label = "3D", linewidths = 1.5)
    ax.set(xlabel = "x", ylabel = "y", zlabel = "z")
    ax.view_init(30, 30)

    fig.suptitle(fname)
    fig.set_tight_layout(True)
    plt.show()


def main():
    for i in range(1, len(sys.argv)):
        plotData(sys.argv[i])


if __name__ == "__main__":
    main()
