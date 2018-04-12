dist_name:=ITK
src_dir:=./src
build_dir:=./build
bin_dir:=./bin
lib_dir:=./lib
share_dir:=./share
install_dir:=..
include_dir:=./include

GIT:=$(shell which git)
TAR:=$(shell which tar)
MAKE:=$(shell which make)
CMAKE:=$(shell which cmake)
CC:=$(shell which cc)
CXX:=$(shell which c++)
AR:=$(shell which ar)
RANLIB:=$(shell which ranlib)

export MAKE CMAKE CC CXX AR RANLIB

checkout_cmd:=$(GIT) clone git://itk.org/ITK.git

cmake_cmd:=$(CMAKE) ../$(src_dir) \
-G "Unix Makefiles" \
-DITK_BUILD_DEFAULT_MODULES=OFF \
-DITKGroup_Core=ON \
-DITKGroup_Filtering=ON \
-DCMAKE_BUILD_TYPE=Release \
-DCMAKE_CXX_COMPILER=$(CXX) \
-DCMAKE_CXX_COMPILER_AR=$(AR) \
-DCMAKE_CXX_COMPILER_RANLIB=$(RANLIB) \
-DCMAKE_CXX_FLAGS="-msse2 -mfpmath=sse" \
-DCMAKE_C_COMPILER=$(CC) \
-DCMAKE_C_COMPILER_AR=$(AR) \
-DCMAKE_C_COMPILER_RANLIB=$(RANLIB) \
-DCMAKE_C_FLAGS="-msse2 -mfpmath=sse" \
-DCMAKE_VERBOSE_MAKEFILE=ON \
-DCMAKE_INSTALL_PREFIX=$(install_dir)

# default target is first target
all: checkout config build test install

.PHONY: nuke checkout config build config test install

nuke:
	rm -rf $(dist_name) $(src_dir) $(build_dir) $(bin_dir) $(lib_dir) $(include_dir) $(share_dir)

checkout:
	$(checkout_cmd)
	mv $(dist_name) $(src_dir)
	# explicitly get latest release, v4.13.0
	(cd $(src_dir) && git checkout v4.13.0)

config:
	mkdir -p $(build_dir)
	(cd $(build_dir) && rm -Rf CMakeCache.txt CMakeFiles/ CMakeTmp/)
	(cd $(build_dir) && $(cmake_cmd))

build:
	(cd $(build_dir) && $(MAKE) clean)
	(cd $(build_dir) && $(MAKE) -j8)

test:
	(cd $(build_dir) && $(MAKE) test)

install:
	(cd $(build_dir) && $(MAKE) install)

debug:
	@echo "checkout command = $(checkout_cmd)" 
	@echo "cmake command = $(cmake_cmd)" 


