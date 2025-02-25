# Rules for compiling the PDF from the LaTeX sources and displaying the output

# Variables
###############################################################################
### Git tags for different revisions
SUBMITTED = submitted-gji
### Documents to build
PDF = paper/preprint.pdf paper/manuscript.pdf paper/manuscript-diff-$(SUBMITTED).pdf paper/misc/cover-letter.pdf
### File Types (for dependencies)
TEX = $(filter-out $(PDF:.pdf=.tex), $(wildcard paper/*.tex))
TEXVARS = $(wildcard paper/variables/*.tex)
BIB = $(wildcard paper/*.bib)
FIG = $(wildcard paper/figures/*)

# Main targets
###############################################################################
preprint: paper/preprint.pdf

manuscript: paper/manuscript.pdf

diff-submitted: paper/manuscript-diff-$(SUBMITTED).pdf

letter: paper/misc/cover-letter.pdf

all: $(PDF)

clean:
	rm -rf $(PDF) paper/manuscript-diff-*.tex paper/_diff

format:
	black code/

lock: conda-lock.yml

# Rules for building PDFs
###############################################################################
%.pdf: %.tex $(TEX) $(BIB) $(FIG)
	tectonic -X compile $<

paper/manuscript-diff-$(SUBMITTED).tex: paper/manuscript.tex $(TEX) $(BIB) $(FIG)
	cd paper \
		&& latexdiff-vc --git --flatten --force --dir=_diff --rev $(SUBMITTED) `basename $<` \
		&& cd .. \
		&& mv paper/_diff/*.tex $@

paper/misc/cover-letter.pdf: paper/misc/cover-letter.tex paper/information.tex
	tectonic -X compile $<

paper/variables.tex: $(TEXVARS)
	cat $^ > $@

# Rules for generating results and other files
###############################################################################
conda-lock.yml: environment.yml
	conda-lock -f $<

# paper/figures/%.png paper/variables/%.tex &: code/%.ipynb code/euler.py
# 	jupyter execute --inplace --kernel_name=python3 $< && touch paper/figures/$*.png paper/variables/$*.tex
#
# paper/figures/real-data-application.png paper/variables/real-data-application.tex &: code/real-data-application.ipynb code/euler.py data/rio-de-janeiro-magnetic.csv
# 	jupyter execute --inplace --kernel_name=python3 $< && touch paper/figures/real-data-application.png paper/variables/real-data-application.tex
#
# data/rio-de-janeiro-magnetic.csv paper/variables/real-data-preparation.tex &: code/real-data-preparation.ipynb data/raw/1038_XYZ.tar.xz
# 	jupyter execute --inplace --kernel_name=python3 $< && touch data/rio-de-janeiro-magnetic.csv paper/variables/real-data-preparation.tex
