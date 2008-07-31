project(tetgen)
include_directories(${tetgen_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(tetgen_SRCS
predicates.cxx tetgen.cxx
)
add_executable(tetgen-test ${tetgen_SRCS})
target_link_libraries(tetgen-test ${FS_LIBS})

add_definitions (-DTETLIBRARY)
add_library (tetgen ${tetgen_SRCS})


