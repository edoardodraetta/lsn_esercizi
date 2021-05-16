CC = g++
CFLAGS = -Wall -O3 --std=c++11
OBJS = random.o metropolis.o

main.exe : $(OBJS)
	$(CC) $(OBJS) -o main.exe

metropolis.o : metropolis.cpp random.h
random.o : random.cpp

clean :
	rm $(OBJS) main.exe
