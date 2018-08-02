#!/usr/bin/env bash

set -x
# set -ex
# set -o pipefail

#
# the gcc5 mac output is riddled with innocuous deprecation warnings -
# we'll have to deal with these until gcc6, but for now, we should filter the output
# so the travis log doesn't get filled up
#

cmake . -DFS_PACKAGES_DIR="packages" -DBUILD_GUIS=OFF

if [[ "$TRAVIS_OS_NAME" == "osx" ]] ; then
  make 2>&1 | grep -v -e '^/var/folders/*' -e '^[[:space:]]*\.section' -e '^[[:space:]]*\^[[:space:]]*~*'
else
  make mris2rgb VERBOSE=1
  cd mris2rgb
  /usr/lib/ccache/g++   -fopenmp  -Wall -Wno-unused-but-set-variable -Wno-unused-result -Wno-unused-local-typedefs -msse2 -mfpmath=sse -g -O3 -m64 -fdata-sections -ffunction-sections -Wl,--gc-sections   -fdata-sections -ffunction-sections -Wl,--gc-sections -Wl,-Map,ld_map.txt -Wl,--no-demangle CMakeFiles/mris2rgb.dir/mris2rgb.c.o  -o mris2rgb ../utils/libutils.a ../packages/glut/3.7/lib/libglut.a /usr/lib/x86_64-linux-gnu/libX11.so /usr/lib/x86_64-linux-gnu/libXmu.so ../rgb/librgb.a ../hipsstubs/libhipsstubs.a ../dicom/libdicom.a ../log/liblog.a ../unix/libunix.a /usr/lib/x86_64-linux-gnu/libz.so ../packages/xml2/2.7.7/lib/libxml2.a ../packages/jpeg/6b/lib/libjpeg.a ../packages/tiff/3.6.1/lib/libtiff.a ../packages/minc/1.5/lib/libvolume_io.a ../packages/minc/1.5/lib/libminc.a ../packages/netcdf/3.6.0/lib/libnetcdf.a ../packages/itk/5.0.0/lib/libITKLabelMap-5.0.a ../packages/itk/5.0.0/lib/libITKSpatialObjects-5.0.a ../packages/itk/5.0.0/lib/libITKPath-5.0.a ../packages/itk/5.0.0/lib/libITKQuadEdgeMesh-5.0.a ../packages/itk/5.0.0/lib/libITKMesh-5.0.a ../packages/itk/5.0.0/lib/libITKOptimizers-5.0.a ../packages/itk/5.0.0/lib/libITKStatistics-5.0.a ../packages/itk/5.0.0/lib/libitkNetlibSlatec-5.0.a ../packages/itk/5.0.0/lib/libITKPolynomials-5.0.a ../packages/itk/5.0.0/lib/libITKBiasCorrection-5.0.a ../packages/itk/5.0.0/lib/libITKIOBMP-5.0.a ../packages/itk/5.0.0/lib/libITKIOGDCM-5.0.a ../packages/itk/5.0.0/lib/libitkgdcmMSFF-5.0.a ../packages/itk/5.0.0/lib/libitkgdcmDICT-5.0.a ../packages/itk/5.0.0/lib/libitkgdcmIOD-5.0.a ../packages/itk/5.0.0/lib/libITKEXPAT-5.0.a ../packages/itk/5.0.0/lib/libitkgdcmDSED-5.0.a ../packages/itk/5.0.0/lib/libitkgdcmCommon-5.0.a ../packages/itk/5.0.0/lib/libitkgdcmjpeg8-5.0.a ../packages/itk/5.0.0/lib/libitkgdcmjpeg12-5.0.a ../packages/itk/5.0.0/lib/libitkgdcmjpeg16-5.0.a ../packages/itk/5.0.0/lib/libitkgdcmopenjp2-5.0.a ../packages/itk/5.0.0/lib/libitkgdcmcharls-5.0.a ../packages/itk/5.0.0/lib/libitkgdcmuuid-5.0.a ../packages/itk/5.0.0/lib/libITKIOGIPL-5.0.a ../packages/itk/5.0.0/lib/libITKIOJPEG-5.0.a ../packages/itk/5.0.0/lib/libITKIOMeta-5.0.a ../packages/itk/5.0.0/lib/libITKMetaIO-5.0.a ../packages/itk/5.0.0/lib/libITKIONIFTI-5.0.a ../packages/itk/5.0.0/lib/libITKTransform-5.0.a ../packages/itk/5.0.0/lib/libITKniftiio-5.0.a ../packages/itk/5.0.0/lib/libITKznz-5.0.a ../packages/itk/5.0.0/lib/libITKIONRRD-5.0.a ../packages/itk/5.0.0/lib/libITKNrrdIO-5.0.a ../packages/itk/5.0.0/lib/libITKIOPNG-5.0.a ../packages/itk/5.0.0/lib/libitkpng-5.0.a ../packages/itk/5.0.0/lib/libITKIOTIFF-5.0.a ../packages/itk/5.0.0/lib/libitktiff-5.0.a ../packages/itk/5.0.0/lib/libitkzlib-5.0.a ../packages/itk/5.0.0/lib/libitkjpeg-5.0.a ../packages/itk/5.0.0/lib/libITKIOVTK-5.0.a ../packages/itk/5.0.0/lib/libITKIOImageBase-5.0.a ../packages/itk/5.0.0/lib/libITKCommon-5.0.a ../packages/itk/5.0.0/lib/libitkdouble-conversion-5.0.a ../packages/itk/5.0.0/lib/libitksys-5.0.a ../packages/itk/5.0.0/lib/libITKVNLInstantiation-5.0.a ../packages/itk/5.0.0/lib/libitkvnl_algo-5.0.a ../packages/itk/5.0.0/lib/libitkvnl-5.0.a ../packages/itk/5.0.0/lib/libitkv3p_netlib-5.0.a ../packages/itk/5.0.0/lib/libitknetlib-5.0.a ../packages/itk/5.0.0/lib/libitkvcl-5.0.a -lm -lpthread -ldl -lm -lcrypt -lrt ../packages/glut/3.7/lib/libglut.a /usr/lib/x86_64-linux-gnu/libX11.so /usr/lib/x86_64-linux-gnu/libXmu.so /usr/lib/x86_64-linux-gnu/libGLU.so /usr/lib/x86_64-linux-gnu/libGL.so
fi