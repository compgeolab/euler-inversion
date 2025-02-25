# Rules for compiling the PDF from the LaTeX sources and displaying the output

### Documents to build
PDF = paper/preprint.pdf paper/manuscript.pdf paper/manuscript-diff-R1.pdf
### File Types (for dependencies)
TEX = $(filter-out $(PDF:.pdf=.tex), $(wildcard paper/*.tex))
TEXVARS = $(wildcard paper/variables/*.tex)
BIB = $(wildcard paper/*.bib)
FIG = $(wildcard paper/figures/*)
### Git tags for different revisions
SUBMITTED = submitted-gji

preprint: paper/preprint.pdf

manuscript: paper/manuscript.pdf

diff-R1: paper/manuscript-diff-R1.pdf

all: $(PDF)

show: $(PDF)
	@for f in $?; do xdg-open $$f; done

clean:
	rm -f $(PDF) paper/misc/cover-letter.pdf paper/content-diff-*.tex paper/abstract-diff-*.tex \
		paper/*.aux paper/*.log paper/*.bbl paper/*.dvi paper/*.blg paper/*.out

format:
	black code/

lock: conda-lock.yml

%.pdf: %.tex $(TEX) $(BIB) $(FIG)
	tectonic -X compile $<

paper/manuscript-diff-R1.pdf: paper/manuscript-diff-R1.tex paper/content-diff-R1.tex paper/abstract-diff-R1.tex $(TEX) $(BIB) $(FIG)
	tectonic -X compile $<

paper/%-diff-R1.tex: paper/%.tex
	cd paper \
		&& latexdiff-vc --git --flatten --force --dir R1 --rev $(SUBMITTED) `basename $<` \
		&& cd .. \
		&& mv paper/R1/*.tex $@

conda-lock.yml: environment.yml
	conda-lock -f $<

letter: paper/misc/cover-letter.pdf

paper/misc/cover-letter.pdf: paper/misc/cover-letter.tex paper/information.tex
	tectonic -X compile $<

paper/variables.tex: $(TEXVARS)
	cat $^ > $@

# paper/figures/%.png paper/variables/%.tex &: code/%.ipynb code/euler.py
# 	jupyter execute --inplace --kernel_name=python3 $< && touch paper/figures/$*.png paper/variables/$*.tex
#
# paper/figures/real-data-application.png paper/variables/real-data-application.tex &: code/real-data-application.ipynb code/euler.py data/rio-de-janeiro-magnetic.csv
# 	jupyter execute --inplace --kernel_name=python3 $< && touch paper/figures/real-data-application.png paper/variables/real-data-application.tex
#
# data/rio-de-janeiro-magnetic.csv paper/variables/real-data-preparation.tex &: code/real-data-preparation.ipynb data/raw/1038_XYZ.tar.xz
# 	jupyter execute --inplace --kernel_name=python3 $< && touch data/rio-de-janeiro-magnetic.csv paper/variables/real-data-preparation.tex
