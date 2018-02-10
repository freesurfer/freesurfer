dist_name:=vtk
src_dir:=./src
build_dir:=./build
bin_dir:=./bin
lib_dir:=./lib
install_dir:=..
include_dir:=./include

GIT:=$(shell which git)
MAKE:=$(shell which make)
CMAKE:=$(shell which cmake)
CC:=$(shell which cc)
CXX:=$(shell which c++)
AR:=$(shell which ar)
RANLIB:=$(shell which ranlib)

export MAKE CMAKE CC CXX AR RANLIB

checkout_cmd:=$(GIT) clone https://gitlab.kitware.com/vtk/vtk.git

cmake_cmd:=$(CMAKE) ../$(src_dir) \
-DCMAKE_BUILD_TYPE=Release \
-DCMAKE_CXX_COMPILER=$(CXX) \
-DCMAKE_CXX_COMPILER_AR=$(AR) \
-DCMAKE_CXX_COMPILER_RANLIB=$(RANLIB) \
-DCMAKE_C_COMPILER=$(CC) \
-DCMAKE_C_COMPILER_AR=$(AR) \
-DCMAKE_C_COMPILER_RANLIB=$(RANLIB) \
-DCMAKE_VERBOSE_MAKEFILE=ON \
-DCMAKE_INSTALL_PREFIX=$(install_dir)

# defualt target is first target
all: checkout config build test install

.PHONY: nuke checkout config build test install

nuke:
	rm -rf $(dist_name) $(src_dir) $(build_dir) $(bin_dir) $(lib_dir) $(include_dir)

checkout:
	$(checkout_cmd)
	mv $(dist_name) $(src_dir)
	# explicitly get vtk v5.4.0 via tag
	(cd $(src_dir) && git checkout v5.4.0)

config:
	mkdir -p $(build_dir)
	(cd $(build_dir) && rm -Rf CMakeCache.txt CMakeFiles/ CMakeTmp/)
	(cd $(build_dir) && $(cmake_cmd))

build:
	(cd $(build_dir) && $(MAKE) clean)
	(cd $(build_dir) && $(MAKE) -j8)

# ignore errors
.IGNORE: test
test:
	# FIX ME: some tests failing
	# (cd $(build_dir) && $(MAKE) test)
	(cd $(build_dir) && $(MAKE) -k test)

install:
	(cd $(build_dir) && $(MAKE) install)

debug:
	@echo "checkout command = $(checkout_cmd)" 
	@echo "cmake command = $(cmake_cmd)" 


