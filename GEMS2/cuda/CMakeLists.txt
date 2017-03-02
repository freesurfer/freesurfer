SET( cudasrcs cudacontroller.cpp )
LIST(APPEND cudasrcs atlasmeshvisitcountercuda.cpp)

CUDA_ADD_LIBRARY( kvlGEMSCUDA ${cudasrcs} )