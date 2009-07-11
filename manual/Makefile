
PDFLATEX=pdflatex
BIBTEX=bibtex

.PHONY : clean

default: all

all:	cl-genomic-manual.pdf

cl-genomic-manual.pdf: cl-genomic-manual.tex
	${PDFLATEX} cl-genomic-manual.tex
	${BIBTEX} cl-genomic-manual
	${PDFLATEX} cl-genomic-manual.tex

clean:
	rm -f cl-genomic-manual.pdf
	rm -f *.aux *.bbl *.blg *.log *.out
