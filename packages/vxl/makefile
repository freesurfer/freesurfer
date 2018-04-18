dist_name:=vxl
src_dir:=./src
bin_dir:=./bin
build_dir:=./build
lib_dir:=./lib
install_dir:=..
include_dir:=./include
share_dir:=./share

GIT:=$(shell which git)
MAKE:=$(shell which make)
CMAKE:=$(shell which cmake)
CC:=$(shell which cc)
CXX:=$(shell which c++)
AR:=$(shell which ar)
RANLIB:=$(shell which ranlib)

export MAKE CMAKE CC CXX AR RANLIB

checkout_cmd:=$(GIT) clone https://github.com/vxl/vxl.git

cmake_cmd:=$(CMAKE) ../$(src_dir) \
-G "Unix Makefiles" \
-DBOXM2_USE_VOLM=OFF \
-DBUILD_BRL=OFF \
-DBUILD_BRL_PYTHON=OFF \
-DBUILD_CONTRIB=OFF \
-DBUILD_CONVERSIONS=OFF \
-DBUILD_CORE_GEOMETRY=OFF \
-DBUILD_CORE_IMAGING=OFF \
-DBUILD_CORE_NUMERICS=ON \
-DBUILD_CORE_SERIALISATION=OFF \
-DBUILD_CORE_UTILITIES=ON \
-DBUILD_CORE_VIDEO=OFF \
-DBUILD_CUL=OFF \
-DBUILD_DOCUMENTATION=OFF \
-DBUILD_EXAMPLES=OFF \
-DBUILD_FOR_VXL_DASHBOARD=OFF \
-DBUILD_GEL=OFF \
-DBUILD_MUL=OFF \
-DBUILD_MUL_TOOLS=OFF \
-DBUILD_NONDEPRECATED_ONLY=OFF \
-DBUILD_OUL=OFF \
-DBUILD_OXL=OFF \
-DBUILD_PRIP=OFF \
-DBUILD_RPL=OFF \
-DBUILD_TBL=OFF \
-DBUILD_TESTING=ON \
-DBUILD_VGUI=OFF \
-DVNL_CONFIG_CHECK_BOUNDS=ON \
-DVNL_CONFIG_ENABLE_SSE2=OFF \
-DVNL_CONFIG_LEGACY_METHODS=ON \
-DVNL_CONFIG_THREAD_SAFE=ON \
-DVXL_LEGACY_ERROR_REPORTING=OFF \
-DWARN_DEPRECATED=OFF \
-DWARN_DEPRECATED_ABORT=OFF \
-DWARN_DEPRECATED_ONCE=OFF \
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

.PHONY: nuke checkout config build test install

nuke:
	rm -rf $(dist_name) $(src_dir) $(build_dir) $(bin_dir) $(lib_dir) $(include_dir) $(share_dir)

checkout:
	$(checkout_cmd)
	mv $(dist_name) $(src_dir)
	# explicitly get vxl v1.17.0 via tag
	# Latest release is >5 years old, and doesn't build using cmake v3.11/gcc v5.4
	# For now just building from most recent commit
	# (cd $(src_dir) && git checkout v1.17.0)

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


