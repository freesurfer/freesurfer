project(mri_sclimbic_seg)

if (FSPYTHON_INSTALL_TREE)
   set(sclimbic_dest "\.\./fspython/models")
   install_pyscript_fspython_tree(mri_sclimbic_seg)
   install_pyscript_fspython_tree(dsurffe)
else()
   set(sclimbic_dest "models")
   install_pyscript(mri_sclimbic_seg)
   install_pyscript(dsurffe)
endif()

install(FILES sclimbic.fsm+ad.t1.nstd00-50.nstd32-50.h5 sclimbic.ctab sclimbic.volstats.csv DESTINATION ${sclimbic_dest})
install(FILES entowm.fsm31.t1.nstd00-30.nstd21-108.h5 entowm.ctab DESTINATION ${sclimbic_dest})
install(FILES mca-dura.both-lh.nstd21.fhs.h5 mca-dura.ctab  DESTINATION ${sclimbic_dest})
install(FILES vsinus.no-sp.m.all.nstd10-070.h5 DESTINATION ${sclimbic_dest})

