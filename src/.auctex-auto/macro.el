(TeX-add-style-hook
 "macro"
 (lambda ()
   (TeX-add-symbols
    '("myacronym" 2)
    '("myglossaryentry" 4))
   (LaTeX-add-index-entries
    "#2@#1"
    "#1"))
 :latex)

