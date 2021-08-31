# Script for plotting the fitness values of the populations of the NSGA-II and NSGA-III algorithms.

import matplotlib.pyplot as plt

# Filenames
f11 = "nsga2_kur_last.txt"
f12 = "nsga2_kur_sols.txt"

f21 = "nsga3_dtlz1_last.txt"
f22 = "nsga3_dtlz1_sols.txt"

# Read and return data points from file
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

# Plot 2D fitness values
def plotData2D(name1, name2, title):
    x1, y1 = readData2D(name1)
    x2, y2 = readData2D(name2)

    fig = plt.figure()
    ax = fig.add_subplot(1, 2, 1)
    ax.scatter(x1, y1, c = "red", marker = "x", label = "Last population", linewidths = 1.0)
    ax.set(xlabel = "f1", ylabel = "f2")
    ax.set_title(title)
    ax.invert_xaxis()
    ax.invert_yaxis()

    ax = fig.add_subplot(1, 2, 2)
    ax.scatter(x2, y2, c = "blue", marker = "x", label = "Optimal solutions", linewidths = 1.0)
    ax.set(xlabel = "f1", ylabel = "f2")
    ax.set_title(title)
    ax.invert_xaxis()
    ax.invert_yaxis()

    fig.legend(loc = "upper right")

    plt.show()
    return

# Plot 3D fitness values
def plotData3D(name1, name2, title):
    x1, y1, z1 = readData3D(name1)
    x2, y2, z2 = readData3D(name2)

    fig = plt.figure()
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    ax.scatter(x1, y1, z1, c = "r", marker = "o", label = "Last population", linewidths = 1.5)
    ax.set(xlabel = "f1", ylabel = "f2", zlabel = "f3")
    ax.set_title(title)
    ax.invert_xaxis()
    ax.invert_yaxis()
    ax.invert_zaxis()

    ax = fig.add_subplot(1, 2, 2, projection="3d")
    ax.scatter(x2, y2, z2, c = "b", marker = "o", label = "Optimal solutions", linewidths = 1.0)
    ax.set(xlabel = "f1", ylabel = "f2", zlabel = "f3")
    ax.set_title(title)
    ax.invert_xaxis()
    ax.invert_yaxis()
    ax.invert_zaxis()

    fig.legend(loc = "upper right")

    plt.show()
    return

def main():
    plotData2D(f11, f12, "KUR fitness values")
    plotData3D(f21, f22, "DTLZ1 fitness values")

if __name__ == "__main__":
    main()