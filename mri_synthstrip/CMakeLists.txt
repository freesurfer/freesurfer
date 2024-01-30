project(synthstrip)

install_pyscript(mri_synthstrip)
install_symlinks(TYPE files DESTINATION models synthstrip.1.pt synthstrip.nocsf.1.pt)
## fails w/o cuda libs
# add_test_script(NAME mri_synthstrip_test SCRIPT test.sh)
