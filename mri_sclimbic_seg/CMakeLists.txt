project(mri_sclimbic_seg)

install_pyscript(mri_sclimbic_seg)
install(FILES sclimbic.h5 sclimbic.ctab sclimbic.volstats.csv DESTINATION models)
