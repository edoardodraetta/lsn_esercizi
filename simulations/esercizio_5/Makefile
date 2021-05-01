CXX = g++
CXXFLAGS = -Wall -O3 --std=c++11
LDFLAGS =

objects =  main.o random.o metropolis.o

main.exe : $(objects)
	$(CXX) $(objects) -o main.exe

main.o : random.h
random.o : random.cpp
metropolis.o : metropolis.cpp

.PHONY : clean
clean :
	rm *.o main.exe
