CXX = g++
CXXFLAGS = -Wall -O3 --std=c++11
LDFLAGS =

objects =  MolDyn.o main.o statistics.o

main.exe : $(objects)
	$(CXX) $(objects) -o main.exe

main.o : MolDyn.h statistics.h
MolDyn.o : MolDyn.cpp


.PHONY : clean
clean :
	rm *.o main.exe
