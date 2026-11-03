project(nnUNet)

# set ref for nnUNet package destination
set(NNUNET_PKG_DST "${FSPYTHON_INSTALL_PREFIX}/python/packages/nnUNet_v2")

# set ref for wrapper to be installed in bin/
set(BIN_PROGRAMS process_directory_no_upsample.sh)

# install all the package dependencies
install(DIRECTORY package/acvl_utils
    package/batchgenerators
    package/batchgeneratorsv2
    package/documentation
    package/dynamic_network_architectures
    package/fft_conv_pytorch
    package/nnunetv2
    DESTINATION ${NNUNET_PKG_DST}
    )

install(FILES package/LICENSE
    package/pyproject.toml
    package/readme.md
    package/setup.py
    DESTINATION ${NNUNET_PKG_DST}
    )

# install the python script in the packages dir for Eugenio's hack
install_symlinks(process_directory_no_upsample.py TYPE programs DESTINATION ${NNUNET_PKG_DST})

# the segmentation daemon the GUI runs from NNUNET_SCRIPT_DIR; it imports nnunetv2 from this dir
install_symlinks(process_directory_no_upsample_daemon.py TYPE programs DESTINATION ${NNUNET_PKG_DST})

# install the wrapper in bin/
install_symlinks(${BIN_PROGRAMS} TYPE programs DESTINATION bin)

