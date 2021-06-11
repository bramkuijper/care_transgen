# makefile to build the IBM
# to check the analytical care_transgen
# model
#
EXE=care_transgen.exe

$(EXE) : individual.hpp individual.cpp care_transgen.cpp
	g++ -Wall -std=c++11 -O3 -o $(EXE) individual.cpp care_transgen.cpp

.PHONY : clean

clean :
	rm -rf $(EXE)
