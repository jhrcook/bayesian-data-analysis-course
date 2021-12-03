.PHONY:

# Need to specify bash in order for conda activate to work.
SHELL=/bin/bash

help:
	@echo "available commands"
	@echo " - help               : information about available commands"
	@echo " - render             : render the distill website"
	@echo " - notebooks          : run Notes and Assignment Rmd"
	@echo " - build              : run all Rmd and render the site"

render:
	Rscript -e "rmarkdown::render_site()" && rm -r docs/data

notebooks:
	Rscript run-all-notebooks.R

build: notebooks render
