.PHONY: all clean install dev-install test
SHELL = /bin/bash -e

all: install

install:
	@which pip > /dev/null
	@pip freeze|grep 'pbinternal2=='>/dev/null \
      && pip uninstall -y pbinternal2 \
      || echo -n ''
	@pip install ./


clean:
	rm -rf build/;\
	find . -name "*.egg-info" | xargs rm -rf;\
	find . -name "*.pyc" | xargs rm -f;\
	find . -name "*.err" | xargs rm -f;\
	find . -name "*.log" | xargs rm -f;\
	rm -rf dist/

test:
	nosetests --nocapture --nologcapture --verbose tests/*.py
	#nosetests --verbose tests/unit/*
	#cram tests/cram/*.t

pip-install:
	@which pip > /dev/null
	@pip freeze|grep 'pbinternal2=='>/dev/null \
      && pip uninstall -y pbinternal2 \
      || true
	@pip install --no-index \
          --install-option="--install-scripts=$(PREFIX)/bin" ./

emit-tool-contracts:
	mkdir -p tool-contracts
	python ./bin/mh_toy.py emit-tool-contracts -o tool-contracts
	python -m pbinternal2.analysis_tools emit-tool-contracts -o tool-contracts
	python -m pbinternal2.pa_tasks emit-tool-contracts -o tool-contracts
	python -m pbinternal2.tasks.eol_qc emit-tool-contracts -o tool-contracts
	python -m pbinternal2.tasks.loading emit-tool-contracts -o tool-contracts

emit-tcs: emit-tool-contracts

auto-pep8:
	find bin -name "*.py" -exec autopep8 -i --ignore=E501,E265,E731,E402,W292 {} \;
	find pbinternal2 -name "*.py" -exec autopep8 -i --ignore=E501,E265,E731,E402,W292 {} \;
