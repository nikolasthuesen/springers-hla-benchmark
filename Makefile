# You can set ?= variables from the command line.
SHELL = /bin/bash
VIRT ?= virt


.PHONY: install

install $(VIRT):
	python3 -m venv $(VIRT)
	$(VIRT)/bin/pip install --upgrade pip setuptools wheel
	$(VIRT)/bin/pip install -r requirements.txt
	$(VIRT)/bin/python setup.py develop

run_benchmark $(VIRT):
	snakemake --use-singularity -k --cores 10 --snakefile snakemake/Snakefile --configfile snakemake/config.yaml