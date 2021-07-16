CC = g++
CFLAGS = -Wall -O3 --std=c++11
OBJS = random.o metropolis.o

Metro.exe : $(OBJS)
	$(CC) $(OBJS) -o Metro.exe

metropolis.o : metropolis.cpp random.h
random.o : random.cpp

clean :
	rm $(OBJS) Metro.exe
