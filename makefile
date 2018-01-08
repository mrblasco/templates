VIEWER := open -a Skim
EDITOR := open -a BBEdit

slides.pdf: slides.Rmd
	Rscript -e "rmarkdown::render('$<')"

view:
	$(VIEWER) slides.pdf
	
edit:
		$(EDITOR) *.Rmd *.yml
	
clean:
	rm slides.pdf