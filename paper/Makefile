# Rules for compiling the PDF from the LaTeX sources and displaying the output

### Documents to build
PDF = preprint.pdf
#PDF = preprint.pdf manuscript.pdf
### File Types (for dependencies)
TEX = $(filter-out $(PDF:.pdf=.tex), $(wildcard *.tex))
BIB = $(wildcard *.bib)
FIG = $(wildcard figures/*)

all: $(PDF)

%.pdf: %.tex $(TEX) $(BIB) $(FIG)
	tectonic -X compile $<

show: $(PDF)
	for fname in $^ ; do \
		python -c "import sys; import webbrowser; webbrowser.open(sys.argv[1])" $$fname ; \
	done

clean:
	rm -f $(PDF)
