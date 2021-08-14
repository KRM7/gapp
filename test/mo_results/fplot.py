# Script for plotting the fitness values of the population of the NSGA-II algorithm.

import matplotlib.pyplot as plt

# Filenames
f11 = "nsga2_kur_sols.txt"
f12 = "nsga3_kur_sols.txt"

f21 = "nsga2_zdt2_sols.txt"
f22 = "nsga3_zdt2_sols.txt"

f31 = "nsga2_zdt3_sols.txt"
f32 = "nsga3_zdt3_sols.txt"

f41 = "nsga2_zdt6_sols.txt"
f42 = "nsga3_zdt6_sols.txt"

f51 = "nsga2_dtlz1_sols.txt"
f52 = "nsga3_dtlz1_sols.txt"

f61 = "nsga2_dtlz2_sols.txt"
f62 = "nsga3_dtlz2_sols.txt"

# Read and return data points from file
def readData(name):
    with open(name) as f:
        lines = f.readlines()
        x = [float(line.split()[0]) for line in lines]
        y = [float(line.split()[1]) for line in lines]
    return x, y

# Plot the fitness values of the NSGA-II and III for the same problem
def plotData(name1, name2, title):
    x1, y1 = readData(name1)
    x2, y2 = readData(name2)

    fig, ax = plt.subplots()
    ax.scatter(x1, y1, c = "red", marker = "x", label = "NSGA-II", linewidths = 1.0)
    ax.scatter(x2, y2, c = "blue", marker = "x", label = "NSGA-III", linewidths = 1.0)
    ax.set(xlabel = "f1", ylabel = "f2")
    ax.set_title(title)
    ax.invert_xaxis()
    ax.invert_yaxis()
    fig.legend(loc = "upper right")

    plt.show()
    return

# Plot the fitness values of the last population
def plotData3D(name1, name2, title):
    with open(name1) as f:
        lines = f.readlines()
        x1 = [float(line.split()[0]) for line in lines]
        y1 = [float(line.split()[1]) for line in lines]
        z1 = [float(line.split()[2]) for line in lines]

    with open(name2) as f:
        lines = f.readlines()
        x2 = [float(line.split()[0]) for line in lines]
        y2 = [float(line.split()[1]) for line in lines]
        z2 = [float(line.split()[2]) for line in lines]

    fig = plt.figure()
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    ax.scatter(x1, y1, z1, c = "r", marker = "o", label = "NSGA-II", linewidths = 1.5)
    ax.set(xlabel = "f1", ylabel = "f2", zlabel = "f3")
    ax.set_title(title)
    ax.invert_xaxis()
    ax.invert_yaxis()
    ax.invert_zaxis()

    ax = fig.add_subplot(1, 2, 2, projection="3d")
    ax.scatter(x2, y2, z2, c = "b", marker = "o", label = "NSGA-III", linewidths = 1.5)
    ax.set(xlabel = "f1", ylabel = "f2", zlabel = "f3")
    ax.set_title(title)
    ax.invert_xaxis()
    ax.invert_yaxis()
    ax.invert_zaxis()

    fig.legend(loc = "upper right")

    plt.show()
    return

def main():
    plotData(f11, f12, "KUR fitness values")
    plotData(f21, f22, "ZDT2 fitness values")
    plotData(f31, f32, "ZDT3 fitness values")
    plotData(f41, f42, "ZDT6 fitness values")

    plotData3D(f51, f52, "DTLZ1")
    plotData3D(f61, f62, "DTLZ2")

if __name__ == "__main__":
    main()