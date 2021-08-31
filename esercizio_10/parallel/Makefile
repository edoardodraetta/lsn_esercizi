# CXX = g++
CXX = mpic++
CXXFLAGS = -Wall -O2 --std=c++11
OBJS = random.o ParallelSearch.o
LDFLAGS = -larmadillo

all : TSP.exe

TSP.exe : $(OBJS)
	$(CXX) $(OBJS) -o TSP.exe $(LDFLAGS)

ParallelSearch.o : ParallelSearch.cpp ParallelSearch.h random.h
	$(CXX) -c ParallelSearch.cpp -o ParallelSearch.o $(CXXFLAGS)
random.o : random.h
	$(CXX) -c random.cpp -o random.o $(CXXFLAGS)
clean :
	rm *.o TSP.exe
