import math

def main():

    C = 128
    R = 3.442619855899
    V = 9.91256303526217e-3

    #C = 256
    #R = 3.6541528853610088
    #V = 0.00492867323399

    x = [0.0] * (C + 1)

    f = math.exp(-0.5 * R * R)

    x[0] = V / f
    x[1] = R
    x[C] = 0.0

    for i in range(2, C):
        x[i] = math.sqrt(-2.0 * math.log(V / x[i - 1] + f))
        f = math.exp(-0.5 * x[i] * x[i])

    for i in range(0, C + 1):
        print(f"{x[i]},")

if __name__ == "__main__":
    main()
