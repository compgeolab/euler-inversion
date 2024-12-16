# Rules for compiling the PDF from the LaTeX sources and displaying the output

### Documents to build
PDF = paper/preprint.pdf paper/manuscript.pdf
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

format:
	black code/

lock: conda-lock.yml

conda-lock.yml: environment.yml
	conda-lock -f $<

paper/variables.tex: $(TEXVARS)
	cat $^ > $@

paper/figures/%.png paper/variables/%.tex &: code/%.ipynb code/euler.py
	jupyter execute --inplace --kernel_name=python3 $< && touch paper/figures/$*.png paper/variables/$*.tex

paper/figures/real-data-application.png paper/variables/real-data-application.tex &: code/real-data-application.ipynb code/euler.py data/rio-de-janeiro-magnetic.csv
	jupyter execute --inplace --kernel_name=python3 $< && touch paper/figures/real-data-application.png paper/variables/real-data-application.tex

data/rio-de-janeiro-magnetic.csv paper/variables/real-data-preparation.tex &: code/real-data-preparation.ipynb data/raw/1038_XYZ.tar.xz
	jupyter execute --inplace --kernel_name=python3 $< && touch data/rio-de-janeiro-magnetic.csv paper/variables/real-data-preparation.tex
