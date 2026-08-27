project(mri_sclimbic_seg)

if (FSPYTHON_INSTALL_TREE)
   set(sclimbic_dest "${FSPYTHON_INSTALL_PREFIX}/models")
   install_modelpredict_pyscript_fspython_tree(sclimbic_seg)
else()
   set(sclimbic_dest "models")
   install_modelpredict_pyscript(sclimbic_seg)
endif()

install(FILES sclimbic.fsm+ad.t1.nstd00-50.nstd32-50.h5 sclimbic.ctab sclimbic.volstats.csv DESTINATION ${sclimbic_dest})
install(FILES entowm.fsm31.t1.nstd00-30.nstd21-108.h5 entowm.ctab DESTINATION ${sclimbic_dest})
install(FILES mca-dura.both-lh.nstd21.fhs.h5 mca-dura.ctab  DESTINATION ${sclimbic_dest})
install(FILES vsinus.no-sp.m.all.nstd10-070.h5 DESTINATION ${sclimbic_dest})

