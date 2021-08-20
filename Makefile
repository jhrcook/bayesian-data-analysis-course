.PHONY: munge download_data docs

# Need to specify bash in order for conda activate to work.
SHELL=/bin/bash

help:
	@echo "available commands"
	@echo " - help               : information about available commands"
	@echo " - render             : render the distill website"

render:
	Rscript -e "rmarkdown::render_site()"
