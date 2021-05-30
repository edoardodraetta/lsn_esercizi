CXX = g++
CXXFLAGS = -Wall -O2 --std=c++11
OBJS = random.o GeneticAlgorithm.o
LDFLAGS = -larmadillo

TSP.exe : $(OBJS)
	$(CXX) $(OBJS) -o TSP.exe $(LDFLAGS)

GeneticAlgorithm.o : GeneticAlgorithm.cpp GeneticAlgorithm.h random.h
	$(CXX) -c GeneticAlgorithm.cpp -o GeneticAlgorithm.o $(CXXFLAGS)
random.o : random.h
	$(CXX) -c random.cpp -o random.o $(CXXFLAGS)
clean :
	rm *.o TSP.exe
