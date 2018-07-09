project(qdec_glmfit)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(qdec_glmfit qdec_glmfit.cxx)
target_link_libraries(qdec_glmfit utils)

install(TARGETS qdec_glmfit DESTINATION bin)
