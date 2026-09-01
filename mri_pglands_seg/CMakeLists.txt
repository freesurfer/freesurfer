project(mri_pglands_seg)
#add_subdirectory(pglands_seg)
# install the bash wrapper for the main entry point,
if (FSPYTHON_INSTALL_TREE)
    set(pglands_dest "\.\./fspython/models/pglands_seg")
    install_modelpredict_pyscript_fspython_tree(pglands_seg)
else()
    set(pglands_dest "models/pglands_seg")
    install_modelpredict_pyscript(pglands_seg)
endif()

# use a soft link for mri_pglands_seg to be included in the fspylegacy as pglands_seg.py
ADD_CUSTOM_TARGET(mri_pglands_seg ALL COMMAND ${CMAKE_COMMAND} -E create_symlink \"${CMAKE_SOURCE_DIR}/mri_pglands_seg/mri_pglands_seg\" \"${CMAKE_SOURCE_DIR}/python/fspylegacy/pglands_seg.py\")
set_property(DIRECTORY ${CMAKE_SOURCE_DIR}/mri_pglands_seg APPEND PROPERTY ADDITIONAL_CLEAN_FILES ${CMAKE_SOURCE_DIR}/python/fspylegacy/pglands_seg.py)

# install the cLUT, template and model
install_symlinks(TYPE FILES pglands_seg/mni152_label_template.mgz pglands_seg/pglands.ctab pglands_seg/pglands_seg.pth DESTINATION ${pglands_dest})

# unpack and install the qa label volumes
install_tarball(pglands_seg_qa_data.tar.gz DESTINATION ${pglands_dest})
