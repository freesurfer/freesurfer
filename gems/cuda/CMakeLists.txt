set(cudasrcs cudacontroller.cpp)

list(APPEND cudasrcs atlasmeshvisitcountercuda.cpp)
list(APPEND cudasrcs visitcountersimplecudaimpl.cu)
list(APPEND cudasrcs visitcountertetrahedralmeshcudaimpl.cu)

list(APPEND cudasrcs atlasmeshalphadrawercuda.cpp)
list(APPEND cudasrcs atlasmeshalphadrawercudaimpl.cu)

cuda_add_library(kvlGEMSCUDA ${cudasrcs})
