CC = g++
CFLAGS = -Wall -O3 --std=c++11
OBJS = random.o Variational_MC.o

Variational_MC.exe : $(OBJS)
	$(CC) $(OBJS) -o Variational_MC.exe

Variational_MC.o : Variational_MC.cpp random.h

random.o : random.cpp

clean :
	rm $(OBJS) Variation_MC.exe
