# FNM lab 1 Makefile
.SUFFIXES: .c

CC = gcc
#CFLAGS =   -Wall -ansi -pedantic -std=c99 -march=native -fomit-frame-pointer -O2
CFLAGS = -g -Wall -ansi -pedantic -std=c99 -march=native
LOADLIBES = -lgsl -lgslcblas -lm

decode: decode.c
	$(CC) $(CFLAGS) $< -o decode $(LOADLIBES)
