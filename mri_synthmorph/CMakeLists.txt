project(mri_synthmorph)

file(GLOB WEIGHTS "synthmorph.*.h5")

if (FSPYTHON_INSTALL_TREE)
   install_pyscript_fspython_tree(mri_synthmorph)
   install_symlinks_fspython_tree(TYPE files DESTINATION models ${WEIGHTS})
   install_directories_fspython_tree(synthmorph DESTINATION python/packages)
else()
   install_pyscript(mri_synthmorph)
   install_symlinks(TYPE files DESTINATION models ${WEIGHTS})
   install_directories(synthmorph DESTINATION python/packages)
endif()

install_configured(fs-synthmorph-reg DESTINATION bin)

## 08/2024 - currently failing on Intel Mac

if(NOT APPLE)
   add_test_script(NAME mri_synthmorph_test_register SCRIPT test_register.sh)
   add_test_script(NAME mri_synthmorph_test_apply SCRIPT test_apply.sh)
endif()
