project(mri_sclimbic_seg)

install_pyscript(mri_sclimbic_seg)
install(FILES sclimbic.fsm+ad.t1.nstd00-50.nstd32-50.h5 sclimbic.ctab sclimbic.volstats.csv DESTINATION models)
