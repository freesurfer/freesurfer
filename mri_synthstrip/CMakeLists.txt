project(synthstrip)

if (FSPYTHON_INSTALL_TREE)
   install_pyscript_fspython_tree(mri_synthstrip)
   install_symlinks_fspython_tree(TYPE files DESTINATION models synthstrip.1.pt synthstrip.nocsf.1.pt)
else()
   install_pyscript(mri_synthstrip)
   install_symlinks(TYPE files DESTINATION models synthstrip.1.pt synthstrip.nocsf.1.pt)
endif()

# cmake ... -DTEST_WITH_CUDA=ON ...
if(TEST_WITH_CUDA)
   add_test_script(NAME mri_synthstrip_test SCRIPT test.sh)
endif()

