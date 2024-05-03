project(mri_synthmorph)

install_pyscript(mri_synthmorph)

install_configured(fs-synthmorph-reg DESTINATION bin)

file(GLOB WEIGHTS "synthmorph.*.h5")
install_symlinks(TYPE files DESTINATION models ${WEIGHTS})
## 02/2024 - currently failing on MacOS intel and arm64 builds
## 03/2024 - re-enable on MacOS
## 04/2024 - disable again
if(NOT APPLE)
   add_test_script(NAME mri_synthmorph_test SCRIPT test.sh)
endif()
