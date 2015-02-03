# FNM lab 1 Makefile
.SUFFIXES: .tex .gp .eps .pdf .c

CC = gcc
# if optimization should be needed...
#CFLAGS = -Wall -ansi -pedantic -fomit-frame-pointer -O2 -march=native -std=c99
# if debug info counts
CFLAGS = -g -Wall -ansi -pedantic -std=c99 -march=native
LOADLIBES = -lgsl -lgslcblas -lm
LC = pdflatex
LATEX_FLAGS = -shell-escape
GNUPLOT_TO_TEX = gnuplot -e $(subst filename,$@,"set terminal epslatex color solid; set output 'filename'")

# When new plot DO:
# 1: Add to PLOTS list
# 2: Create dependencies list and compilation rule
PLOTS = filtered_with_bandwidth_256 different_filters different_bandwidths

rap: rap.tex $(PLOTS:=.tex)
	$(LC) $(LATEX_FLAGS) rap.tex

# Plots dependencies and compilation rules
filtered_with_bandwidth_256.tex: filtered_with_bandwidth_256.gp amplot.data filtered_256.data
	$(GNUPLOT_TO_TEX) $< # $< means "the first dependency". Here: filtered_with_bandwidth_256.dp

different_filters.tex: different_filters.gp fft_before_filtering.data filtered_fft_domain_64.data Scale_fct_64.data filtered_fft_domain_128.data Scale_fct_128.data filtered_fft_domain_256.data Scale_fct_256.data
	$(GNUPLOT_TO_TEX) $<

different_bandwidths.tex: different_bandwidths.gp filtered_64.data filtered_128.data filtered_256.data filtered_512.data
	$(GNUPLOT_TO_TEX) $<

# _POSIX_C_SOURCE >= 2 needed to get popen from stdio.h
lab1: lab1.c
	$(CC) -D_POSIX_C_SOURCE=2 $(CFLAGS) $< -o lab1 $(LOADLIBES)

decode: decode.c
	$(CC) $(CFLAGS) $< -o decode $(LOADLIBES)

warmup0: warmup0.c
	$(CC) -D_POSIX_C_SOURCE=2 $(CFLAGS) $< -o warmup0 $(LOADLIBES)

warmup1: warmup1.c
	$(CC) -D_POSIX_C_SOURCE=2 $(CFLAGS) $< -o warmup1 $(LOADLIBES)

.PHONY: clean veryclear show

show:
	mupdf rap.pdf

# remove temp files
clean:
	rm -vf $(PLOTS:=.eps) $(PLOTS:=.tex) $(PLOTS:=-eps-converted-to.pdf) *.aux *.toc *.log 

# also remove all compiled results
veryclean:
	rm -vf $(PLOTS:=.eps) $(PLOTS:=.tex) $(PLOTS:=-eps-converted-to.pdf) *.aux *.toc *.log rap.pdf lab1
