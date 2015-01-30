# FNM lab 1 Makefile
.SUFFIXES: .tex .gp .eps .pdf .c

CC = gcc
# if optimization should be needed...
#CFLAGS = -Wall -ansi -pedantic -fomit-frame-pointer -O2 -march=native -std=c99
# if debug info counts
CFLAGS = -g -Wall -ansi -pedantic -std=c99 -march=native
LOADLIBES = -lgsl -lgslcblas -lm



LC = pdflatex
FLAGS = -shell-escape
PLOTS = 
PLOTINCLUDES = $(PLOTS:=.tex)
GNUPLOTSCRIPTS = $(PLOTS:=.gp)
EPSS = $(PLOTS:=.eps)

rap: rap.tex $(PLOTINCLUDES)
	$(LC) $(FLAGS) rap.tex

# Always make eps and tex when gnuplotting through Make
%.tex: %.gp
	gnuplot -e $(subst filename,$@,"set terminal epslatex color solid; set output 'filename'") $<

# _POSIX_C_SOURCE >= 2 needed to get popen from stdio.h
lab1: lab1.c
	$(CC) -D_POSIX_C_SOURCE=2 $(CFLAGS) $< -o lab1 $(LOADLIBES)

warmup0: warmup0.c
	$(CC) -D_POSIX_C_SOURCE=2 $(CFLAGS) $< -o warmup0 $(LOADLIBES)


.PHONY: clean veryclear show

show:
	mupdf rap.pdf

# remove temp files
clean:
	rm -vf $(EPSS) $(PLOTINCLUDES) *.aux *.toc *.log $(PLOTS:=-eps-converted-to.pdf)

# also remove all compiled results
veryclean:
	rm -vf $(EPSS) $(PLOTINCLUDES) *.aux *.toc *.log $(PLOTS:=-eps-converted-to.pdf) rap.pdf lab1
