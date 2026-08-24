rm *.out *.aux *.blg *.log

pdflatex ms.tex
bibtex ms
pdflatex ms.tex
pdflatex ms.tex

rm *.out *.aux *.blg *.log

#cp ms.pdf ../STATIC/.

cp ms.pdf STATIC/.
