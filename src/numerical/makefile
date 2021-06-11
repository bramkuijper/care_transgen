# make file to compile the care transgen numerical iterations
EXE=care_transgen.exe
CXX=g++
CXXFLAGS=-Wall -O3 -std=c++17 -ggdb
CPP=care_transgen.cpp calculate.cpp
HPP=care_transgen.hpp

$(EXE) : $(CPP) $(HPP)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(CPP)


clean :
	rm -rf $(EXE)
