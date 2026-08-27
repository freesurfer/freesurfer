project(mri_claustrum_seg)

# handle FSPYTHON tree split install
if (FSPYTHON_INSTALL_TREE)
    set(claustrum_python "${FSPYTHON_INSTALL_PREFIX}/python")
    set(claustrum_models "${FSPYTHON_INSTALL_PREFIX}/models")
else()
    set(claustrum_python "python")
    set(claustrum_models "models")
endif()

# bash wrapper always installed the same way, fspython agnostic
install(PROGRAMS mri_claustrum_seg DESTINATION bin)

# install the model file
if (FSPYTHON_INSTALL_TREE)
   install_modelpredict_pyscript_fspython_tree(seg_claustrum)
   install_symlinks_fspython_tree(TYPE files DESTINATION "${claustrum_models}" model/claustrum_seg_20250616.h5)
 else()
   install_modelpredict_pyscript(seg_claustrum)   
   install_symlinks(TYPE files DESTINATION "${claustrum_models}" model/claustrum_seg_20250616.h5)
endif()

## isntall the atlas dir
install_tarball(claustrum_seg.tar.gz DESTINATION average)
