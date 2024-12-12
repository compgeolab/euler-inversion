# Rules for compiling the PDF from the LaTeX sources and displaying the output

### Documents to build
PDF = paper/preprint.pdf
#PDF = preprint.pdf manuscript.pdf
### File Types (for dependencies)
TEX = $(filter-out $(PDF:.pdf=.tex), $(wildcard paper/*.tex))
TEXVARS = $(wildcard paper/variables/*.tex)
BIB = $(wildcard paper/*.bib)
FIG = $(wildcard paper/figures/*)

all: $(PDF)

%.pdf: %.tex $(TEX) $(BIB) $(FIG)
	tectonic -X compile $<

show: $(PDF)
	xdg-open $<

clean:
	rm -f $(PDF)

lock:
	# Create lock files for the current version of the environment
	conda-lock -f environment.yml -p osx-64 -p osx-arm64 -p linux-64 -p win-64
	# Use this instead
	conda-lock render -p linux-64
	conda create -n my-locked-env --file conda-linux-64.lock


paper/variables.tex: $(TEXVARS)
	cat $^ > $@

paper/figures/%.png: code/%.ipynb code/euler.py
	jupyter execute --inplace --kernel_name=python3 $<
	# Because jupyter execute modifies the notebook last
	touch $@
	echo ""

data/rio-de-janeiro-magnetic.csv: code/real-data-preparation.ipynb data/1038_XYZ.tar.xz
	jupyter execute --inplace --kernel_name=python3 $<
	# Because jupyter execute modifies the notebook last
	touch $@

format:
	black code/
