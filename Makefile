APPNAME='tdbsumstat'
TARGETS=build clean dependencies deploy install test uninstall
VERSION=$(shell grep version pyproject.toml|cut -f2 -d'"')


all:
	@echo "Try one of: ${TARGETS}"

build_conda_lock_files:
	conda-lock -k explicit -f base_environment.yml --conda mamba -p 'linux-64' -p 'osx-arm64'

build: clean dependencies
	poetry build

clean:
	find . -name '*.pyc' -delete
	find . -type d -name '__pycache__' -exec rm -rf {} +
	rm -rf dist build

dependencies:
	poetry install --no-root

deploy:
	poetry install

install: build
	pip install dist/*.whl --force-reinstall

editable_install: build
	pip install --editable .

tag:
	git tag v${VERSION}

uninstall:
	pip uninstall -y ${APPNAME}
