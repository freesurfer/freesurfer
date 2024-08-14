project(mri_synthmorph)

file(GLOB WEIGHTS "synthmorph.*.h5")

if (FSPYTHON_INSTALL_TREE)
   install_pyscript_fspython_tree(mri_synthmorph)
   install_symlinks_fspython_tree(TYPE files DESTINATION models ${WEIGHTS})
else()
   install_pyscript(mri_synthmorph)
   install_symlinks(TYPE files DESTINATION models ${WEIGHTS})
endif()

install_configured(fs-synthmorph-reg DESTINATION bin)

## 08/2024 - currently failing on Intel Mac

if(NOT APPLE)
   add_test_script(NAME mri_synthmorph_test SCRIPT test.sh)
elseif(APPLE_ARM64)
   add_test_script(NAME mri_synthmorph_test SCRIPT test.sh)
endif()

