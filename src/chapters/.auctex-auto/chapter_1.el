(TeX-add-style-hook
 "chapter_1"
 (lambda ()
   (LaTeX-add-labels
    "chapter:introduction"
    "eq:1")
   (LaTeX-add-index-entries
    "inline"))
 :latex)

