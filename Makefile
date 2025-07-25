# settings

PY := python3
ENV_NAME := venv

# computed

ACTIV_ENV := . $(ENV_NAME)/bin/activate
IN_ENV := $(ACTIV_ENV) &&

# includes
# -include *.mk

py_init:
	@echo "\n[[ creating python environment if absent ]]"
	$(PY) -m venv $(ENV_NAME)

py_dep: py_init
	@echo "\n[[ installing/upgrading python dependencies ]]"
	$(IN_ENV) python3 -m pip install --upgrade pip
	$(IN_ENV) pip install --upgrade -r requirements.txt

setup: py_init py_dep

rm_cache:
	@echo "\n[[ removing all pycache folders ]]"
	find . -name __pycache__ -prune -print -exec rm -rf {} \;

rm_temp:
	@echo "\n[[ removing generated temporary files ]]"
	rm -f yosys_graph.log

rm_pyenv:
	@echo "\n[[ removing the virtual python environment ]]"
	rm -rf $(ENV_NAME)

clean: rm_pyenv rm_cache rm_temp