CXX = g++
CXXFLAGS = -Wall -O2 --std=c++11
OBJS = random.o anneal.o
LDFLAGS = -larmadillo

all : TSP.exe

TSP.exe : $(OBJS)
	$(CXX) $(OBJS) -o TSP.exe $(LDFLAGS)

anneal.o : anneal.cpp anneal.h random.h
	$(CXX) -c anneal.cpp -o anneal.o $(CXXFLAGS)
random.o : random.h
	$(CXX) -c random.cpp -o random.o $(CXXFLAGS)
clean :
	rm *.o TSP.exe
