(TeX-add-style-hook
 "preamble"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("natbib" "numbers" "sort&compress") ("hyperref" "colorlinks=true" "urlcolor=blue" "anchorcolor=blue" "citecolor=blue" "filecolor=blue" "linkcolor=blue" "menucolor=blue" "pagecolor=blue" "linktocpage=true" "pdfproducer=medialab" "pdfa=true") ("glossaries" "toc" "nonumberlist" "nopostdot") ("geometry" "outer=36.02988mm" "inner=25.47697mm" "top=36.0317mm" "bottom=50.9565mm") ("caption" "labelfont={bf,singlespacing}" "textfont={singlespacing}" "justification={justified,RaggedRight}" "singlelinecheck=false" "margin=0pt" "figurewithin=chapter" "tablewithin=chapter") ("appendix" "titletoc") ("quotchap" "helvetica") ("background" "contents={}") ("datetime" "yyyymmdd" "hhmmss")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "amsmath"
    "epsfig"
    "graphicx"
    "natbib"
    "color"
    "hyperref"
    "url"
    "makeidx"
    "glossaries"
    "geometry"
    "lipsum"
    "caption"
    "appendix"
    "quotchap"
    "background"
    "datetime")
   (TeX-add-symbols
    "markblankpages"
    "symmetricmargin"
    "bookmargin"
    "DraftText"
    "draft"
    "archivalpapernote"))
 :latex)

