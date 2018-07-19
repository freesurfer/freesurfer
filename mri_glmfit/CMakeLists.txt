project(mri_glmfit)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_glmfit mri_glmfit.c)
target_link_libraries(mri_glmfit utils fsgdf)

add_test_script(NAME mri_glmfit_test SCRIPT test_mri_glmfit DEPENDS mri_glmfit)

install_wrapped(TARGETS mri_glmfit DESTINATION bin)
