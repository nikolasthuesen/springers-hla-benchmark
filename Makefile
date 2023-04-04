# You can set ?= variables from the command line.
SHELL = /bin/bash
VIRT ?= virt


install:
	python3 -m venv $(VIRT)
	$(VIRT)/bin/pip install --upgrade pip setuptools wheel
	$(VIRT)/bin/pip install -r requirements.txt
	$(VIRT)/bin/pip install --editable .
	$(VIRT)/bin/python setup.py develop


run_benchmark:
	snakemake \
		--use-singularity \
		--singularity-args "-B $(PWD):$(PWD)" \
		-k \
		--verbose \
		--latency-wait 20 \
		--cores 16 \
		--snakefile snakemake/Snakefile \
		--configfile snakemake/config.yaml 
	
	
run_slim_benchmark:
	snakemake \
		--use-singularity \
		--singularity-args "-B $(PWD):$(PWD)" \
		-k \
		--verbose \
		--latency-wait 20 \
		--cores 24 \
		--snakefile snakemake/Snakefile_slim \
		--configfile snakemake/config.yaml

clear_snakemake_cache:
	rm -r .snakemake/auxiliary/ .snakemake/conda* .snakemake/incomplete/ .snakemake/lo* .snakemake/metadata/ .snakemake/shadow/