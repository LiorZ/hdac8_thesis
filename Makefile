all: thesis.pdf
	cp thesis.pdf ~/Dropbox/Thesis_Lior/.
thesis.pdf: thesis.tex thesis.aux
	pdflatex thesis.tex
thesis.tex: thesis.rst
	rst2latex --latex-preamble "\usepackage[hmargin=3cm,vmargin=3.5cm]{geometry} \usepackage{setspace}" thesis.rst thesis.tex
	 
