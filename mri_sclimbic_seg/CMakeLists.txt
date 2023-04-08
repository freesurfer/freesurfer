project(mri_sclimbic_seg)

install_pyscript(mri_sclimbic_seg)
install(FILES sclimbic.fsm+ad.t1.nstd00-50.nstd32-50.h5 sclimbic.ctab sclimbic.volstats.csv DESTINATION models)
install(FILES entowm.fsm31.t1.nstd00-30.nstd21-108.h5 entowm.ctab DESTINATION models)
