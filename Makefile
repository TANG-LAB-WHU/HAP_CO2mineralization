PYTHON ?= python3

.PHONY: setup check test evidence-local-check paper paper-pdf preview clean

setup:
	bash scripts/bootstrap_macos.sh --install

check:
	$(PYTHON) scripts/check_manuscript.py
	$(PYTHON) -m unittest discover -s tests -p 'test_*.py'

test:
	$(PYTHON) -m unittest discover -s tests -p 'test_*.py'

evidence-local-check:
	$(PYTHON) scripts/check_manuscript.py --check-local-pdfs

paper: check
	quarto render paper --to html

paper-pdf: check
	quarto render paper --to pdf

preview: check
	quarto preview paper --to html

clean:
	rm -rf paper/_output paper/.quarto
