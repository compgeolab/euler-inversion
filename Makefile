# Rules for compiling the PDF from the LaTeX sources and displaying the output

### Documents to build
PDF = paper/preprint.pdf
#PDF = preprint.pdf manuscript.pdf
### File Types (for dependencies)
TEX = $(filter-out $(PDF:.pdf=.tex), $(wildcard paper/*.tex))
BIB = $(wildcard paper/*.bib)
FIG = $(wildcard paper/figures/*)

all: $(PDF)

%.pdf: %.tex $(TEX) $(BIB) $(FIG)
	tectonic -X compile $<

show: $(PDF)
	xdg-open $<

clean:
	rm -f $(PDF)
