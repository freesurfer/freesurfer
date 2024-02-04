project(synthstrip)

install_pyscript(mri_synthstrip)
install_symlinks(TYPE files DESTINATION models synthstrip.1.pt synthstrip.nocsf.1.pt)
# cmake ... -DTEST_WITH_CUDA=ON ...
if(TEST_WITH_CUDA)
   add_test_script(NAME mri_synthstrip_test SCRIPT test.sh)
endif()

