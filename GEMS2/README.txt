I used cmake-3.5.2, which has a nice Matlab Mex compilation option without any further hassles. 

In all CMake commands below, make sure to set the CMAKE_BUILD_TYPE to "Release" (optimized and thus faster code)

In general, if you're interested in the Matlab Mex files then it's easiest to compile everything using shared libraries everywhere. 
Otherwise, on amd64 (x86_64) platform all ITK etc. libraries will need to be built with Position Independent Code, i.e., the -fPIC 
switch for gcc, because we will be building a shared library (mex-files really are renamed shared libraries). In order to do that, 
set the following compilation flags in CMake: CMAKE_CXX_FLAGS=-fPIC and CMAKE_C_FLAGS=-fPIC.

At the Martinos center, the Matlab set-up is very convoluted and CMake doesn't it find everything it needs automatically. The trick is to set the Matlab path straight to it's actual location:

  /usr/pubsw/common/matlab/8.4/
  
Other than that, I used the following versions:

* InsightToolkit-4.9.1 (shared libs on, V3 compatibility on)
  (InsightToolkit-3.20.0.tar.gz doesn't compile with my current compiler)
  V3 compatibility is needed for kvlRegisterer.h/cxx since you nowadays have to call "Update()" to start registration, and itkAffine3DTransform.h/cxx (also only needed for kvlRegisterer and kvlResampler) needs nowadays different Jacobian interfaces. Nothing major - we never use Jacobian so could just throw an exception if someone attempts to use - but I simply don't have time for it now. In any case, do we really still need kvlRegisterer and kvlResampler in the long run?
 
CMake will need the ITK_DIR environment variable set to the directory containing ITKConfig.cmake

* gmm-5.0 (GMM++ C++ templated library for sparse matrix solvers; used only in Levenberg-Marquardt. If nobody is using this, maybe we can just throw it out in the future?)

  ./configure 
  make

Note that this is a header library, so to actually compile anything you will need to run
  make check

* tetgen1.4.2 (used in one of the initial stages of the atlas building, so as to start with a sparse mesh in the background thereby saving mesh simplification time)

    tar xvfz tetgen1.4.2.tar.gz 
    cd tetgen1.4.2

    Usually you'd now type "make tetlib", but if you've switched shared libs on, you need to manually do it with the -fPIC option specified:

    g++ -O0 -c -fPIC predicates.cxx
    g++ -g -Wall -DSELF_CHECK -DTETLIBRARY -fPIC -c tetgen.cxx
    ar r libtet.a tetgen.o predicates.o

For GUI:   
   
  * VTK5.10.1 (with lots of warnings using new CMake, but seems to work fine)
  Using older version because VTK7 and even VTK6 has changed pipeline architecture from "b->SetInput( a->GetOutput() )" into "b->SetInputConnection( a->GetOutputPort() )". That's of course trivial to change in kvlImageViewer.cxx (which is the only VTK-using class), but there are also non-trial ITK-to-VTK pipeline connectors for which I really just don't have the time to adjust right now.

  * fltk-1.1.10 (I also tried the newest version fltk-1.3.3 because it supports CMake builds, but sadly (1) had some link problems for me when trying to use CMake, and (2) even without CMake I got same very strange link problems in my own GUI code)

     ./configure --enable-shared
     make -j 8
