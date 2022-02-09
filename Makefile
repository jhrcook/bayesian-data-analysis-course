.PHONY:

# Need to specify bash in order for conda activate to work.
SHELL=/bin/bash

help:
	@echo "available commands"
	@echo " - help               : information about available commands"
	@echo " - book              : run all notebooks and build the gitbook"

book:
	Rscript -e 'bookdown::render_book("index.Rmd", "bookdown::gitbook")'
