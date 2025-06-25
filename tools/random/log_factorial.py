import math

def approx_log_factorial(k):
    log_sqrt_2pi = 0.9189385332046727;

    k1 = k + 1.0;
    k1_inv = 1.0 / (k + 1.0);
    k1_inv_sq = k1_inv * k1_inv;

    return log_sqrt_2pi + (k + 0.5) * math.log(k1) - k1 + (1.0 / 12.0 - (1.0 / 360.0 - (1.0 / 1260.0 * k1_inv_sq)) * k1_inv_sq) * k1_inv;


def main():
    for i in range(0, 32):
        exact = math.log(math.factorial(i))
        approx = approx_log_factorial(i)
        print(f"k={i}, value: {exact}, approx: {approx}, error: {math.fabs(exact-approx)}, ulp: {math.ulp(exact)}")


if __name__ == "__main__":
    main()
