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

os:=$(shell uname -s | tr -s '[A-Z]' '[a-z]')
kern:=$(shell uname -v | tr -s '[A-Z]' '[a-z]')

UBUNTU:=false
CENTOS:=false
ifneq (,$(findstring ubuntu,$(kern)))
   UBUNTU:=true
endif
ifneq (,$(wildcard /etc/redhat-release))
   CENTOS:=true
endif

# Need at least cmake version 3.3 and gcc 5.5
#
# CentOS7 has 4.X gcc and devtoolset-4 has gcc 5.3 which does not work either.
# There is no devtoolset-5 with a gcc > 5.3, so install devtoolset-6 with gcc6.3
# yum install cmake3.x86_64
# yum install centos-release-scl
# yum-config-manager --enable rhel-server-rhscl-7-rpms
# yum install devtoolset-6
#
# on Mac OS with homebrew, "brew install gcc@5" which is gcc 5.5
#
# on Ubuntu 16.04, nothing to add as already have gcc 5.4 and cmake 3.5.1

ifeq ($(os),linux)
   ifeq ($(CENTOS),true)
      cmake_bin_dir:=/usr/bin
      CMAKE:=$(cmake_bin_dir)/cmake3
      gcc_bin_dir=/opt/rh/devtoolset-6/root/usr/bin
      CC:=$(gcc_bin_dir)/gcc
      CXX:=$(gcc_bin_dir)/g++
      PATH:=$(cmake_bin_dir):$(gcc_bin_dir):$(PATH)
   else ifeq ($(UBUNTU),true)
      cmake_bin_dir:=/usr/bin
      CMAKE:=$(cmake_bin_dir)/cmake
      gcc_bin_dir=/usr/bin
      CC:=$(gcc_bin_dir)/gcc
      CXX:=$(gcc_bin_dir)/g++
   endif
else ifeq ($(os),darwin)
   cmake_bin_dir:=/usr/local/bin
   CMAKE:=$(cmake_bin_dir)/cmake
   gcc_bin_dir:=/usr/local/bin
   CC:=$(gcc_bin_dir)/gcc-5
   CXX:=$(gcc_bin_dir)/g++-5
   PATH:=$(cmake_bin_dir):$(gcc_bin_dir):$(PATH)
endif

# vxl git source most recent commit hash as of 4/18/2018 (no recent tag will build)
checkout_tag:=5d00d94afade

AR:=$(shell which ar)
RANLIB:=$(shell which ranlib)

export PATH MAKE CMAKE CC CXX AR RANLIB

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
# Don't build anything if the compiler cannot be found

ifeq (,$(wildcard $(CC)))

# susbtitue "error" for "warning" below if want build to fail with non-zero status if vxl cannot build
.PHONY: all
all:
	$(warning Cannot stat required compiler $(CC) for $(os) vxl build, skipping vxl and continuing build.)

else

.PHONY: nuke checkout config build test install
all: checkout config build test install

nuke:
	rm -rf $(dist_name) $(src_dir) $(build_dir) $(bin_dir) $(lib_dir) $(include_dir) $(share_dir)

checkout:
	$(checkout_cmd)
	mv $(dist_name) $(src_dir)
	(cd $(src_dir) && git checkout $(checkout_tag))

config:
	mkdir -p $(build_dir)
	(cd $(build_dir) && rm -Rf CMakeCache.txt CMakeFiles/ CMakeTmp/)
	(cd $(build_dir) && $(cmake_cmd))

build:
	(cd $(build_dir) && $(MAKE) clean)
	(cd $(build_dir) && $(MAKE) -j4)

test:
	(cd $(build_dir) && $(MAKE) test)

install:
	(cd $(build_dir) && $(MAKE) install)

debug:
	@echo "checkout command = $(checkout_cmd)" 
	@echo "cmake command = $(cmake_cmd)" 

endif

