#!/usr/bin/env bash

# build Tcl/Tk 8.4.6, Tix 8.1.4 and BLT 2.4 in such a way
# that they all play nicely with each other, and with the freesurfer
# tools such as tkmedit, tksurfer and tkregister2

set -e

export CC=$(which gcc)
export CXX=$(which g++)

if [ "$#" != "1" ] ; then echo "error: usage: build.sh <prefix>" && exit 1 ; fi
INSTALL_DIR="$1"

# COMBO_DIR must be set to where installation is desired
export COMBO_DIR="${INSTALL_DIR}"

# bug fix for problem that prevents tksurfer from globbing directories on 64 bit NFS file systems
[ "$(uname -s)" == "Linux" ] && export CFLAGS="-DHAVE_STRUCT_DIRENT64=1"

export MACOSX_DEPLOYMENT_TARGET="10.4"
export PATH="${COMBO_DIR}/src/tcl8.4.6/unix:${COMBO_DIR}/src/tk8.4.6/unix:${PATH}"
export LD_LIBRARY_PATH="${COMBO_DIR}/src/tcl8.4.6/unix"
export ENABLE_SYMBOLS=""


# -- build TCL and TK --

# tcl
cd tcl8.4.6/unix
./configure ${ENABLE_SYMBOLS} \
    --prefix=${COMBO_DIR} \
    --exec-prefix=${COMBO_DIR} \
    --mandir='${prefix}/share/man'
make -j8
make install
if [ -e "libtcl8.4g.so" ] ; then
    rm -f libtcl8.4.so
    ln -s libtcl8.4g.so libtcl8.4.so
fi

# tk
cd ../../tk8.4.6/unix
./configure ${ENABLE_SYMBOLS} \
    --prefix=${COMBO_DIR} \
    --exec-prefix=${COMBO_DIR} \
    --mandir='${prefix}/share/man'
make -j8
if [ -e "libtk8.4g.so" ] ; then
    rm -f libtk8.4.so
    ln -s libtk8.4g.so libtk8.4.so
fi

# install
cd ../../tcl8.4.6/unix
make genstubs
make install
cd ../../tk8.4.6/unix
make install

[ -e "${COMBO_DIR}/lib/libtcl.so.0"   ] && rm -f               ${COMBO_DIR}/lib/libtcl.so.0
[ -e "${COMBO_DIR}/lib/libtcl8.4g.so" ] && rm -f               ${COMBO_DIR}/lib/libtcl8.4.so
[ -e "${COMBO_DIR}/lib/libtcl8.4g.so" ] && ln -s libtcl8.4g.so ${COMBO_DIR}/lib/libtcl8.4.so
[ -e "${COMBO_DIR}/lib/libtcl8.4.so"  ] && ln -s libtcl8.4.so  ${COMBO_DIR}/lib/libtcl.so.0
[ -e "${COMBO_DIR}/lib/libtk.so.0"    ] && rm -f               ${COMBO_DIR}/lib/libtk.so.0
[ -e "${COMBO_DIR}/lib/libtk8.4g.so"  ] && rm -f               ${COMBO_DIR}/lib/libtk8.4.so
[ -e "${COMBO_DIR}/lib/libtk8.4g.so"  ] && ln -s libtk8.4g.so  ${COMBO_DIR}/lib/libtk8.4.so
[ -e "${COMBO_DIR}/lib/libtk8.4.so"   ] && ln -s libtk8.4.so   ${COMBO_DIR}/lib/libtk.so.0

# These are needed for FindTCL to find this built instance rather than /usr/bin etc.
[ -e "${COMBO_DIR}/bin/tclsh8.4"      ] && ln -s tclsh8.4      ${COMBO_DIR}/bin/tclsh
[ -e "${COMBO_DIR}/bin/wish8.4"       ] && ln -s wish8.4       ${COMBO_DIR}/bin/wish
[ -e "${COMBO_DIR}/lib/libtcl8.4.so"  ] && ln -s libtcl8.4.so  ${COMBO_DIR}/lib/libtcl.so
[ -e "${COMBO_DIR}/lib/libtk8.4.so"   ] && ln -s libtk8.4.so   ${COMBO_DIR}/lib/libtk.so

unset CFLAGS


# -- build TIX --

cd ../../tix-8.1.4/unix

./configure --prefix="${COMBO_DIR}" \
    --with-tclconfig="${COMBO_DIR}/lib" \
    --with-tkconfig="${COMBO_DIR}/lib" \
    --with-tclinclude=${COMBO_DIR}/include \
    --with-tkinclude=${COMBO_DIR}/include

cd tk8.4

export TIXLIB_BUILD_TYPE="enable-shared"
# on the Mac, just build static module, not the shared, as the shared has a
# compatibility-verison of 0.0.0, which causes trouble
[ "$(uname -s)" == "Linux" ] && export TIXLIB_BUILD_TYPE="enable-sam"

./configure --prefix="${COMBO_DIR}" \
    --$TIXLIB_BUILD_TYPE \
    --with-tclconfig="${COMBO_DIR}/lib" \
    --with-tkconfig="${COMBO_DIR}/lib" \
    --with-tclinclude=${COMBO_DIR}/include \
    --with-tkinclude=${COMBO_DIR}/include

make -j8 all TCL_SRC_DIR=../../../tcl8.4.6 TK_SRC_DIR=../../../tk8.4.6
make -j8 prefix=${COMBO_DIR} \
    MAN_DIR=${COMBO_DIR}/share/man install \
    TCL_SRC_DIR=../../../tcl8.4.6 \
    TK_SRC_DIR=../../../tk8.4.6 \
    RUN_TCLSH=${COMBO_DIR}/bin/tclsh8.4

/usr/bin/install -d ${COMBO_DIR}/share/doc/tix

cd ..
make install

[ -e "${COMBO_DIR}/lib/libtix.so.0"      ] && rm -f ${COMBO_DIR}/lib/libtix.so.0
[ -e "${COMBO_DIR}/lib/libtix8.1.8.4.so" ] && ln -s libtix8.1.8.4.so ${COMBO_DIR}/lib/libtix.so.0
[ -e "${COMBO_DIR}/lib/libtix8.1.8.3.so" ] && rm -f ${COMBO_DIR}/lib/libtix8.1.8.3.so
[ -e "${COMBO_DIR}/lib/libtix8.1.8.4.so" ] && ln -s libtix8.1.8.4.so ${COMBO_DIR}/lib/libtix8.1.8.3.so
[ -e "${COMBO_DIR}/lib/libtix8.1.8.4.a"  ] && chmod 755 ${COMBO_DIR}/lib/libtix*.a
[ -e "${COMBO_DIR}/lib/libtix8.1.8.4.a"  ] && ranlib ${COMBO_DIR}/lib/libtix*.a


# -- build BLT --

cd ../../blt2.4z

[ "$(uname -s)" == "Darwin" ] && export USRLIBDIR="/usr/X11R6/lib"

./configure --prefix=${COMBO_DIR} \
  --mandir=${COMBO_DIR}/share/man \
  --enable-shared \
  --with-tcl=${COMBO_DIR} \
  --with-tk=${COMBO_DIR} \
  --enable-jpeg=${COMBO_DIR} \
  --x-includes=/usr/X11R6/include \
  --x-libraries=/usr/X11R6/lib \
  --with-cflags="-fno-common"

make -j8
make install prefix=${COMBO_DIR}

[ -e "${COMBO_DIR}/lib/libBLT.2.dylib"   ] && rm -f ${COMBO_DIR}/lib/libBLT.2.dylib
[ -e "${COMBO_DIR}/lib/libBLT.2.4.dylib" ] && ln -s libBLT.2.4.dylib ${COMBO_DIR}/lib/libBLT.2.dylib
[ -e "${COMBO_DIR}/lib/libBLT24.a"       ] && chmod 755 ${COMBO_DIR}/lib/libBLT24.a
[ -e "${COMBO_DIR}/lib/libBLT24.a"       ] && ranlib ${COMBO_DIR}/lib/libBLT24.a
