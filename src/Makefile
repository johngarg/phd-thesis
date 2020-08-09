ROOT_DIRECTORY := $(shell pwd)
DIRS := $(sort $(dir $(wildcard $(ROOT_DIRECTORY)/*/$(ROOT_DIRECTORY)/*/*/ )))
CLEAN_DIRS := $(sort $(dir $(shell find . -name '*.tex')))
TIMESTAMP := $(shell date +'%Y_%m_%d-%H_%M_%S')
ROOT_TEX := $(wildcard *.tex)
CHAPTERS := $(basename $(notdir $(wildcard chapters/*.tex)))
CHAPTER_RELEASE := $(addsuffix _release,$(CHAPTERS))

BIB = bibtex
GLOSSARIES = makeglossaries
INDICES = makeindex
COMPILER = 	xelatex -pdf -pdflatex="pdflatex -interaction=nonstopmode" -use-make $<

clean:
	@for d in $(CLEAN_DIRS); do\
		cd $$d; pwd; latexmk -c -silent; cd $(ROOT_DIRECTORY);\
	done
	@rm -vf chapter_*.tex

cleanall: clean
	@rm -vf *.pdf
	@find . -name '*.bbl' | xargs rm -v
	@rm -f *.glsdefs

all: main

chapters: $(CHAPTERS)

main: main.tex
	cp $@.pdf $@_view.pdf
	@rm $@.pdf
	cp $@_view.pdf ../thesis.pdf

release: main
	cp $<_view.pdf ../drafts/all_$(TIMESTAMP).pdf

$(ROOT_TEX):
	rm -f *.glsdefs
	xelatex -pdf -pdflatex="pdflatex -interaction=nonstopmode" $@
	$(BIB) $(basename $@)
	$(GLOSSARIES) $(basename $@)
	$(INDICES) $(basename $@)
	xelatex -pdf -pdflatex="pdflatex -interaction=nonstopmode" $@
	xelatex -pdf -pdflatex="pdflatex -interaction=nonstopmode" $@

$(CHAPTERS): %: chapters/%.tex
	./skeleton chapters/$@.tex

$(CHAPTER_RELEASE): %_release: %
	cp $<_view.pdf ../drafts/$<_$(TIMESTAMP).pdf

.PHONY: clean cleanall all main chapters release $(ROOT_TEX)

#---
# catch all 
#---
help: 
	@bash -c "echo Available targets:; \
            make -pn | grep -e '^[a-z0-9-]\+:' | sed 's/:.*//g' | \
            sed 's/\(.*\)/  \1/g' | sort"

# run "make help" if no target is supplied
.DEFAULT_GOAL := help

# a catch-all message for any targets that do not exist
%::
	@$(ECHO) Target \`$@\' not found. Use \"make help\" to list all targets.
