(TeX-add-style-hook "writeup"
 (lambda ()
    (LaTeX-add-environments
     "lemma")
    (TeX-add-symbols
     '("hwsolution" 1)
     '("hwproblem" 2)
     "hwheading"
     "class"
     "date"
     "prof"
     "name"
     "title"
     "R"
     "F"
     "N"
     "Q"
     "Z"
     "C"
     "a"
     "b"
     "d"
     "e"
     "w"
     "Span"
     "char"
     "im"
     "Hom"
     "deg")
    (TeX-run-style-hooks
     "fancyhdr"
     "setspace"
     "tikz"
     "wasysym"
     "listings"
     "color"
     "url"
     "amsthm"
     "amsfonts"
     "amssymb"
     "amsmath"
     "graphicx"
     "enumerate"
     "geometry"
     "latex2e"
     "art12"
     "article"
     "12pt")))

