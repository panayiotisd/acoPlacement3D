# Makefile for acoPlacement
# Copyright (C) 2015  Panayiotis Danassis (panos_dan@hotmail.com)
CC = gcc
OPTIM_FLAGS = -O3
WARN_FLAGS = -Wall
LDLIBS = -lm

all: CFLAGS = $(WARN_FLAGS) $(OPTIM_FLAGS) -c
all: acoPlacement

debug: CFLAGS = $(WARN_FLAGS) -c
debug: CC += -DDEBUG -g
debug: acoPlacement

print: CFLAGS = $(WARN_FLAGS) $(OPTIM_FLAGS) -c
print: CC += -DPRINT
print: acoPlacement

parallel: CFLAGS = $(WARN_FLAGS) $(OPTIM_FLAGS) -c
parallel: CC += -fopenmp -DPARALLEL -DPRINT
parallel: acoPlacement

acoPlacement: acoPlacement.o ants.o inOut.o utilities.o timer.o
	$(CC) acoPlacement.o ants.o inOut.o utilities.o timer.o $(LDLIBS) -o acoPlacement3D

acoPlacement.o: acoPlacement.c
	$(CC) $(CFLAGS) acoPlacement.c

ants.o: ants.c
	$(CC) $(CFLAGS) ants.c

inOut.o: inOut.c
	$(CC) $(CFLAGS) inOut.c

utilities.o: utilities.c
	$(CC) $(CFLAGS) utilities.c

timer.o: timer.c
	$(CC) $(CFLAGS) timer.c

clean:
	rm -rf *.o acoPlacement3D