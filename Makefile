SOURCES = $(wildcard *.c)
OBJS = $(SOURCES:.c=.o)

CC = gcc
LD = gcc

CFLAGS = -O3 -std=c99 -lm -march=native -ftree-vectorize


all: 
	$(LD) $(SOURCES) $(CFLAGS) -o test

clean: 
	rm -rf *.o test
