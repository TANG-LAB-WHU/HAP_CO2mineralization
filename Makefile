PYTHON ?= python3

.PHONY: setup check paper paper-pdf preview clean

setup:
	bash scripts/bootstrap_macos.sh --install

check:
	$(PYTHON) scripts/check_manuscript.py

paper: check
	quarto render paper --to html

paper-pdf: check
	quarto render paper --to pdf

preview: check
	quarto preview paper --to html

clean:
	rm -rf paper/_output paper/.quarto

