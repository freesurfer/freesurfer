dist_name:=tetgen1.5.0
src_dir:=./src
build_dir:=./build
bin_dir:=./bin
lib_dir:=./lib
install_dir:=..

GIT:=$(shell which git)
TAR:=$(shell which tar)
MAKE:=$(shell which make)
CMAKE:=$(shell which cmake)
CC:=$(shell which cc)
CXX:=$(shell which c++)
AR:=$(shell which ar)
RANLIB:=$(shell which ranlib)

export MAKE CMAKE CC CXX AR RANLIB

checkout_cmd:=$(TAR) zxpf tetgen1.5.0.tar.gz

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
	rm -rf $(dist_name) $(src_dir) $(build_dir) $(bin_dir) $(lib_dir)

checkout:
	$(checkout_cmd)
	mv $(dist_name) $(src_dir)

config:
	mkdir -p $(build_dir)
	(cd $(build_dir) && rm -Rf CMakeCache.txt CMakeFiles/ CMakeTmp/)
	(cd $(build_dir) && $(cmake_cmd))

build:
	(cd $(build_dir) && $(MAKE) clean)
	(cd $(build_dir) && $(MAKE) -j8)

test:
	# (cd $(build_dir) && $(MAKE) -k test)

install:
	# (cd $(build_dir) && $(MAKE) install)
	mkdir -p ./bin ./lib
	cp -p -f ./build/tetgen ./bin/
	cp -p -f ./build/libtet.a ./lib/

debug:
	@echo "checkout command = $(checkout_cmd)" 
	@echo "cmake command = $(cmake_cmd)" 


