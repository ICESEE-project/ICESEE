# Dual ISSM build with external package setup and ninja support

ISSM_CLONE_URL := https://github.com/ISSMteam/ISSM.git
ISSM_SRC_DIR := $(PWD)/ISSM
MATLAB_BUILD := $(PWD)/ISSM-matlab
PYTHON_BUILD := $(PWD)/ISSM-python
EXTERNAL_DIR := $(ISSM_SRC_DIR)/externalpackages
PYTHON_VENV := $(HOME)/venv-firedrake
PYTHON_SITE := $(PYTHON_VENV)/lib/python3.11/site-packages/numpy
PYTHON_BASE := /opt/homebrew/opt/python@3.11/Frameworks/Python.framework/Versions/3.11
MATLAB_PATH := /Applications/MATLAB_R2024a.app

.PHONY: all setup external matlab python clean distclean

all: setup external matlab python

setup:
	@if [ ! -d "$(ISSM_SRC_DIR)" ]; then \
		echo "Cloning ISSM..."; \
		git clone $(ISSM_CLONE_URL) $(ISSM_SRC_DIR); \
	fi
	@echo "Adding ISSM environment to ~/.zshrc..."
	@grep -q 'export ISSM_DIR=' ~/.zshrc || echo 'export ISSM_DIR=$(ISSM_SRC_DIR)' >> ~/.zshrc
	@grep -q 'source $$ISSM_DIR/etc/environment.sh' ~/.zshrc || echo 'source $$ISSM_DIR/etc/environment.sh' >> ~/.zshrc

external:
	cd $(EXTERNAL_DIR)/autotools && ./install-mac.sh
	cd $(EXTERNAL_DIR)/cmake && ./install.sh
	cd $(EXTERNAL_DIR)/petsc && ./install-3.22-mac.sh
	cd $(EXTERNAL_DIR)/triangle && ./install-mac.sh
	cd $(EXTERNAL_DIR)/m1qn3 && ./install-mac.sh
	source $(ISSM_SRC_DIR)/etc/environment.sh
	python3 -m pip install --user numpy scipy matplotlib nose
	cd $(ISSM_SRC_DIR) && autoreconf -ivf

matlab:
	@echo "Building ISSM with MATLAB interface..."
	rm -rf $(MATLAB_BUILD)
	mkdir -p $(MATLAB_BUILD)/execution
	rsync -a --exclude='execution' --exclude='.git' $(ISSM_SRC_DIR)/ $(MATLAB_BUILD)/
	cd $(MATLAB_BUILD) && autoreconf -ivf
	cd $(MATLAB_BUILD) && ./configure \
		--prefix=$(MATLAB_BUILD)/install \
		--with-numthreads=2 \
		--with-matlab-dir=$(MATLAB_PATH) \
		--with-fortran-lib="-L/opt/homebrew/Cellar/gcc/14.2.0_1/lib/gcc/14 -lgfortran" \
		--with-mpi-include=$(EXTERNAL_DIR)/petsc/install/include \
		--with-mpi-libflags="-L$(EXTERNAL_DIR)/petsc/install/lib -lmpi -lmpicxx -lmpifort" \
		--with-metis-dir=$(EXTERNAL_DIR)/petsc/install \
		--with-parmetis-dir=$(EXTERNAL_DIR)/petsc/install \
		--with-blas-lapack-dir=$(EXTERNAL_DIR)/petsc/install \
		--with-scalapack-dir=$(EXTERNAL_DIR)/petsc/install \
		--with-mumps-dir=$(EXTERNAL_DIR)/petsc/install \
		--with-petsc-dir=$(EXTERNAL_DIR)/petsc/install \
		--with-triangle-dir=$(EXTERNAL_DIR)/triangle/install \
		--with-m1qn3-dir=$(EXTERNAL_DIR)/m1qn3/install \
		
	cd $(MATLAB_BUILD) && make -j$(shell sysctl -n hw.ncpu) install

python:
	@echo "Building ISSM with Python interface ..."
	rm -rf $(PYTHON_BUILD)
	mkdir -p $(PYTHON_BUILD)/execution
	rsync -a --exclude='execution' --exclude='.git' $(ISSM_SRC_DIR)/ $(PYTHON_BUILD)/
	cd $(PYTHON_BUILD) && autoreconf -ivf
	cd $(PYTHON_BUILD) && ./configure \
		--prefix=$(PYTHON_BUILD)/install \
		--with-python-version=3.11 \
		--with-python-dir=$(PYTHON_BASE) \
		--with-python-numpy-dir=$(PYTHON_SITE) \
		--with-fortran-lib="-L/opt/homebrew/Cellar/gcc/14.2.0_1/lib/gcc/14 -lgfortran" \
		--with-mpi-include=$(EXTERNAL_DIR)/petsc/install/include \
		--with-mpi-libflags="-L$(EXTERNAL_DIR)/petsc/install/lib -lmpi -lmpicxx -lmpifort" \
		--with-metis-dir=$(EXTERNAL_DIR)/petsc/install \
		--with-parmetis-dir=$(EXTERNAL_DIR)/petsc/install \
		--with-blas-lapack-dir=$(EXTERNAL_DIR)/petsc/install \
		--with-scalapack-dir=$(EXTERNAL_DIR)/petsc/install \
		--with-mumps-dir=$(EXTERNAL_DIR)/petsc/install \
		--with-petsc-dir=$(EXTERNAL_DIR)/petsc/install \
		--with-triangle-dir=$(EXTERNAL_DIR)/triangle/install \
		--with-m1qn3-dir=$(EXTERNAL_DIR)/m1qn3/install \

	cd $(PYTHON_BUILD) && make -j$(shell sysctl -n hw.ncpu) install

clean:
	rm -rf $(MATLAB_BUILD) $(PYTHON_BUILD)

distclean: clean
	rm -rf $(EXTERNAL_DIR)/petsc/install \
	       $(EXTERNAL_DIR)/triangle/install \
	       $(EXTERNAL_DIR)/m1qn3/install

.PHONY: test-matlab test-python

test-matlab:
	@echo "Running ISSM MATLAB test suite..."
	@/bin/bash -c '\
		export ISSM_DIR=$(MATLAB_BUILD) && \
		source $(ISSM_DIR)/etc/environment.sh && \
		matlab -r "addpath ${ISSM_DIR}/bin ${ISSM_DIR}/lib" '		

test-python:
	@echo "Running ISSM Python interface check..."
	@/bin/bash -c '\
		source ~/venv-firedrake/bin/activate && \
		export ISSM_DIR=$(PYTHON_BUILD) && \
		python3 -c "\
import os, sys; \
sys.path.append(os.path.join(os.getenv(\"ISSM_DIR\"), \"bin\")); \
sys.path.append(os.path.join(os.getenv(\"ISSM_DIR\"), \"lib\")); \
from issmversion import issmversion; "\
	'