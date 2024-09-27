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

paper/figures/%.png: code/%.ipynb code/euler.py
	jupyter execute --inplace --kernel_name=python3 $<
	# Because jupyter execute modifies the notebook last
	touch $@
	echo ""

data/rio-de-janeiro-magnetic.nc: code/real-data-prepare-grid.ipynb data/1038_XYZ.tar.xz
	jupyter execute --inplace --kernel_name=python3 $<
	# Because jupyter execute modifies the notebook last
	touch $@
	echo ""
