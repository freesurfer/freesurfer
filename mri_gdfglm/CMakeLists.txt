project(mri_gdfglm)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_gdfglm mri_gdfglm.cpp)
target_link_libraries(mri_gdfglm utils)

install(TARGETS mri_gdfglm DESTINATION bin)
