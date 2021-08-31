CXX := g++
CC := g++

WARN_FLAGS := -Wall -Wextra -pedantic
DEBUG_FLAGS := -std=c++20 -Og -ggdb $(WARN_FLAGS)
FLAGS := -std=c++20 -DNDEBUG -O2 -fwhole-program $(WARN_FLAGS)

SRC := test/benchmark.cpp
BIN := benchmark


debug: $(SRC)
	$(CXX) $(DEBUG_FLAGS) $(SRC) -o $(BIN)

benchmark: $(SRC)
	$(CXX) $(FLAGS) $(SRC) -o $(BIN)

ex1 : examples/1_intro_binary_so.cpp
	$(CXX) $(FLAGS) examples/1_intro_binary_so.cpp -o bin/binary_so

ex2 : examples/2_rcga_so.cpp
	$(CXX) $(FLAGS) examples/2_rcga_so.cpp -o bin/rcga_so

ex3 : examples/3_permutations_so.cpp
	$(CXX) $(FLAGS) examples/3_permutations_so.cpp -o bin/permutation_so

ex4 : examples/4_integer_so.cpp
	$(CXX) $(FLAGS) examples/4_integer_so.cpp -o bin/integer_so

ex5: examples/5_nsga2.cpp
	$(CXX) $(FLAGS) examples/5_nsga2.cpp -o bin/nsga2

ex6: examples/6_nsga3.cpp
	$(CXX) $(FLAGS) examples/6_nsga3.cpp -o bin/nsga3

ex7: examples/7_custom_operators.cpp
	$(CXX) $(FLAGS) examples/7_custom_operators.cpp -o bin/custom_operators

clean:
	$(RM) bin/*.exe
	$(RM) *.exe