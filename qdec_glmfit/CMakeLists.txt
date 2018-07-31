project(qdec_glmfit)

include_directories(${FS_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/qdecproject)

add_executable(qdec_glmfit qdec_glmfit.cxx)
target_link_libraries(qdec_glmfit qdecproject)
install(TARGETS qdec_glmfit DESTINATION bin)
