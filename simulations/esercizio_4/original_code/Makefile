
CXX = g++
CXXFLAGS = -Wall -O3 --std=c++11
LDFLAGS =

objects =  MolDyn_NVE.o

MolDyn_NVE.exe : $(objects)
	$(CXX) $(objects) -o MolDyn_NVE.exe

MolDyn_NVE.o : MolDyn_NVE.cpp

.PHONY : clean
clean :
	rm *.o MolDyn_NVE.exe

