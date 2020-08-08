(TeX-add-style-hook
 "thesis"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "unicode-math"
    "enumerate"
    "setspace"
    "fancyhdr")
   (TeX-add-symbols
    '("psection" 1)
    '("summary" 1)
    '("submissionyear" 1)
    '("submissionmonth" 1)
    '("university" 1)
    '("department" 1)
    "smallcapstitle"
    "frontmatterheadings"
    "mainmatterheadings"
    "journalpaperlist"
    "conferencepaperlist"
    "patentlist"
    "patentlistsingleitem"
    "makedeclaration"
    "footnotesize"
    "footnoterule"
    "cleardoublepage")
   (LaTeX-add-environments
    "dedication"
    "abstract"
    "preface"
    "publications"
    "acknowledgements"
    "synopsis")
   (LaTeX-add-lengths
    "bindmargin"
    "othermargin"))
 :latex)

