COMPILER = xelatex -pdf -pdflatex="pdflatex -interactive=nonstopmode" -use-make $<
BIB = bibtex
GLOSSARIES = makeglossaries
INDICES = makeindex
PROJNAME = main
RM = rm -f *.nav *.snm *.fls *.log *.fdb_latexmk *.aux *.synctex.gz *.vrb *.bbl *.glsdefs
RMDIR = rm -rf *.auctex-auto


.PHONY: $(PROJNAME).pdf all clean

all: $(PROJNAME).pdf

$(PROJNAME).pdf: $(PROJNAME).tex
	$(COMPILER) $(PROJNAME).tex
	$(BIB) $(PROJNAME)
	$(GLOSSARIES) $(PROJNAME)
	$(INDICES) $(PROJNAME)
	$(COMPILER) $(PROJNAME).tex
	$(COMPILER) $(PROJNAME).tex

clean:
	latexmk -c

cleanall:
	latexmk -C; $(RM); $(RMDIR)
