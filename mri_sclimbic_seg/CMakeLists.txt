project(mri_sclimbic_seg)

if (FSPYTHON_INSTALL_TREE)
   set(sclimbic_dest "${FSPYTHON_INSTALL_PREFIX}/models")
   install_modelpredict_pyscript_fspython_tree(sclimbic_seg)
else()
   set(sclimbic_dest "models")
   install_modelpredict_pyscript(sclimbic_seg)
endif()

# use a soft link for mri_sclimbic_seg to be included in the fspylegacy as sclimbic_seg.py
ADD_CUSTOM_TARGET(mri_sclimbic_seg ALL COMMAND ${CMAKE_COMMAND} -E create_symlink \"${CMAKE_SOURCE_DIR}/mri_sclimbic_seg/mri_sclimbic_seg\" \"${CMAKE_SOURCE_DIR}/python/fspylegacy/sclimbic_seg.py\")
set_property(DIRECTORY ${CMAKE_SOURCE_DIR}/mri_sclimbic_seg APPEND PROPERTY ADDITIONAL_CLEAN_FILES ${CMAKE_SOURCE_DIR}/python/fspylegacy/sclimbic_seg.py)

install(FILES sclimbic.fsm+ad.t1.nstd00-50.nstd32-50.h5 sclimbic.ctab sclimbic.volstats.csv DESTINATION ${sclimbic_dest})
install(FILES entowm.fsm31.t1.nstd00-30.nstd21-108.h5 entowm.ctab DESTINATION ${sclimbic_dest})
install(FILES mca-dura.both-lh.nstd21.fhs.h5 mca-dura.ctab  DESTINATION ${sclimbic_dest})
install(FILES vsinus.no-sp.m.all.nstd10-070.h5 DESTINATION ${sclimbic_dest})

