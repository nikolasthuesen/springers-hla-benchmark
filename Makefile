# You can set ?= variables from the command line.
SHELL = /bin/bash
VIRT ?= virt


.PHONY: install

install $(VIRT):
	python3 -m venv $(VIRT)
	$(VIRT)/bin/pip install --upgrade pip setuptools wheel
	$(VIRT)/bin/pip install -r requirements.txt
	$(VIRT)/bin/pip install --editable .
	$(VIRT)/bin/python setup.py develop


run_benchmark $(VIRT):
	snakemake \
		--use-singularity \
		--singularity-args "-B $(PWD):$(PWD)" \
		-k \
		--verbose \
		--latency-wait 20 \
		--cores 16 \
		--snakefile snakemake/Snakefile \
		--configfile snakemake/config.yaml 
	
run_slim_benchmark $(VIRT):
	snakemake \
		--use-singularity \
		--singularity-args "-B $(PWD):$(PWD)" \
		-k \
		--verbose \
		--latency-wait 20 \
		--cores 24 \
		--snakefile snakemake/Snakefile_slim \
		--configfile snakemake/config.yaml