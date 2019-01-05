project(mri_glmfit)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_glmfit mri_glmfit.c)
target_link_libraries(mri_glmfit utils fsgdf)

add_test_script(NAME mri_glmfit_test SCRIPT test_mri_glmfit DEPENDS mri_glmfit)
configure_file(teststats.tar.gz . COPYONLY)

install(TARGETS mri_glmfit DESTINATION bin RENAME mri_glmfit)
