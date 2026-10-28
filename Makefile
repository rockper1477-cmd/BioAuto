CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Isrc

all: bin/bioauto bin/tests

bin/bioauto: src/main.cpp $(wildcard src/*.hpp)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ src/main.cpp

bin/tests: tests/tests.cpp $(wildcard src/*.hpp)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ tests/tests.cpp

test: bin/tests
	./bin/tests

clean:
	rm -rf bin
.PHONY: all test clean
