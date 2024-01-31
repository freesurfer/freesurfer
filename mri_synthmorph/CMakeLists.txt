project(mri_synthmorph)

install_pyscript(mri_synthmorph)

file(GLOB WEIGHTS "synthmorph.*.h5")
install_symlinks(TYPE files DESTINATION models ${WEIGHTS})
add_test_script(NAME mri_synthmorph_test SCRIPT test.sh)
