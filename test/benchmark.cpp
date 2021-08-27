#include "time_test.h"

#include "binary_tests.h"
#include "real_tests.h"
#include "permutational_tests.h"
#include "integer_tests.h"

#include "nsga2_tests.h"
#include "nsga3_tests.h"

#include <cstdio>

int main()
{
    binaryRastriginTest();
    binaryRosenbrockTest();
    binarySchwefelTest();
    binaryGriewankTest();
    binaryAckleyTest();

    realRastriginTest();
    realRosenbrockTest();
    realSchwefelTest();
    realGriewankTest();
    realAckleyTest();

    perm52Test();
    perm124Test();
    perm226Test();
    perm439Test();

    integerTest1();
    integerTest2();

    nsga2KurTest();
    nsga2Zdt2Test();
    nsga2Zdt3Test();
    nsga2Zdt6Test();

    nsga3KurTest();
    nsga3Zdt2Test();
    nsga3Zdt3Test();
    nsga3Zdt6Test();

    nsga2Dtlz1Test();
    nsga2Dtlz2Test();

    nsga3Dtlz1Test();
    nsga3Dtlz2Test();

    timeGA();

    std::getchar();

    return 0;
}