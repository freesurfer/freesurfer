project(mri_babyseg)

if (FSPYTHON_INSTALL_TREE)
    install_pyscript_fspython_tree(mri_babyseg)
else()
    install_pyscript(mri_babyseg)
endif()

# install babyseg files
install_symlinks(babyseg/babyseg.v1.pt TYPE files DESTINATION average/babyseg)
install(FILES babyseg/baby.loss.json DESTINATION average/babyseg)
install(FILES babyseg/babyseg.v1.json DESTINATION average/babyseg)
install(FILES babyseg/defaults.json DESTINATION average/babyseg)
install(FILES babyseg/BabySegColorLUT.txt DESTINATION average/babyseg)