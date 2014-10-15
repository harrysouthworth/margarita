all:
	pdflatex  margarita.tex
	bibtex margarita
	pdflatex  margarita.tex
	pdflatex margarita.tex
	pdflatex safety.tex
	bibtex safety
	pdflatex safety.tex
	pdflatex safety.tex
	rm *.log *.aux *.bbl *.blg *.tex
	mv *.pdf ../inst/doc
