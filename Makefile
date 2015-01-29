.SUFFIXES: .tex .gp .eps .pdf
LC = pdflatex
FLAGS = -shell-escape
PLOTS = 
PLOTINCLUDES = $(PLOTS:=.tex)
GNUPLOTSCRIPTS = $(PLOTS:=.gp)
EPSS = $(PLOTS:=.eps)

rap: rap.tex $(PLOTINCLUDES) framsida.pdf
	$(LC) $(FLAGS) rap.tex

# Always make eps and tex when gnuplotting through Make
%.tex: %.gp
	gnuplot -e $(subst filename,$@,"set terminal epslatex color solid; set output 'filename'") $<

framsida.pdf: framsida.tex
	$(LC) $(FLAGS) framsida.tex

.PHONY: clean show

show:
	mupdf rap.pdf

clean:
	rm -vf $(EPSS) $(PLOTINCLUDES) *.aux *.toc *.log framsida.pdf $(PLOTS:=-eps-converted-to.pdf)
