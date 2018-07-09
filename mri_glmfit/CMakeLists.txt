project(mri_glmfit)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_glmfit mri_glmfit.c)
target_link_libraries(mri_glmfit utils)

install(TARGETS mri_glmfit DESTINATION bin)
