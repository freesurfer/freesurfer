project(synthstrip)

# use a soft link for synthstrip to be included in the fspylegacy as synthstrip.py
ADD_CUSTOM_TARGET(synthstrip ALL COMMAND ${CMAKE_COMMAND} -E create_symlink \"${CMAKE_SOURCE_DIR}/mri_synthstrip/mri_synthstrip\" \"${CMAKE_SOURCE_DIR}/python/fspylegacy/synthstrip.py\")
set_property(DIRECTORY ${CMAKE_SOURCE_DIR}/mri_synthstrip APPEND PROPERTY ADDITIONAL_CLEAN_FILES ${CMAKE_SOURCE_DIR}/python/fspylegacy/synthstrip.py)

if (FSPYTHON_INSTALL_TREE)
   install_modelpredict_pyscript_fspython_tree(synthstrip)
   install_symlinks_fspython_tree(TYPE files DESTINATION models synthstrip.1.pt synthstrip.nocsf.1.pt)
else()
   install_modelpredict_pyscript(synthstrip)
   install_symlinks(TYPE files DESTINATION models synthstrip.1.pt synthstrip.nocsf.1.pt)
endif()

# cmake ... -DTEST_WITH_CUDA=ON ...
if(TEST_WITH_CUDA)
   add_test_script(NAME mri_synthstrip_test SCRIPT test.sh)
endif()

