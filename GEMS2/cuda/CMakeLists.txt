SET( cudasrcs cudacontroller.cpp )
LIST(APPEND cudasrcs atlasmeshvisitcountercuda.cpp)
LIST(APPEND cudasrcs visitcountersimplecudaimpl.cu)

CUDA_ADD_LIBRARY( kvlGEMSCUDA ${cudasrcs} )