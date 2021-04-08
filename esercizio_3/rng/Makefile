CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o GeometricBrownianMotion.o statistics.o
	$(CC) random.o statistics.o GeometricBrownianMotion.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
statistics.o : statistics.cpp statistics.h
	$(CC) -c statistics.cpp -o statistics.o $(CFLAGS)
GeometricBrownianMotion.o : GeometricBrownianMotion.cpp GeometricBrownianMotion.h
	$(CC) -c GeometricBrownianMotion.cpp -o GeometricBrownianMotion.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
