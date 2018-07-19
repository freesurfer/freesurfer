SET( cudasrcs cudacontroller.cpp )

LIST(APPEND cudasrcs atlasmeshvisitcountercuda.cpp)
LIST(APPEND cudasrcs visitcountersimplecudaimpl.cu)
LIST(APPEND cudasrcs visitcountertetrahedralmeshcudaimpl.cu)

LIST(APPEND cudasrcs atlasmeshalphadrawercuda.cpp)
LIST(APPEND cudasrcs atlasmeshalphadrawercudaimpl.cu)

CUDA_ADD_LIBRARY( kvlGEMSCUDA ${cudasrcs} )