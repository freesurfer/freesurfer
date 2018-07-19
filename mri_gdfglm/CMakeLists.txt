project(mri_gdfglm)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_gdfglm mri_gdfglm.c)
target_link_libraries(mri_gdfglm utils fsgdf)

install(TARGETS mri_gdfglm DESTINATION bin)
