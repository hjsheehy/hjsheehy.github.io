compile:
	bash ../shell/main.sh
	pdflatex main.tex; \
	bibtex main.tex; \
	pdflatex main.tex; \
	pdflatex main.tex; \
	bash ../shell/rename.sh
	rm main.tex
	rm -f \
		main.aux \
		main.acn \
		main.acr \
		main.alg \
		main.glg \
		main.glo \
		main.gls \
		main.ist \
		main.lof \
		main.lot \
		main.nlg \
		main.not \
		main.ntn \
		main.bbl \
		main.blg \
		main.log \
		main.toc \
		main.out \
		main.glsdefs;
