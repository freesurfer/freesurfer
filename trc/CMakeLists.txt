project(trc)

if(PETSC_FOUND)

  include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${FS_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/fem_elastic
    SYSTEM
    ${CMAKE_SOURCE_DIR}/packages/tetgen
    ${ITK_INCLUDE_DIRS}
    ${PETSC_INCLUDE_DIR}
  )

  # dmri_spline
  add_executable(dmri_spline dmri_spline.cxx spline.cxx)
  target_link_libraries(dmri_spline utils)
  install(TARGETS dmri_spline DESTINATION bin)

  # dmri_forrest
  add_executable(dmri_forrest dmri_forrest.cxx forrest.cxx)
  target_link_libraries(dmri_forrest utils)
  install(TARGETS dmri_forrest DESTINATION bin)

  # dmri_mergepaths
  add_executable(dmri_mergepaths dmri_mergepaths.cxx)
  target_link_libraries(dmri_mergepaths utils)
  install(TARGETS dmri_mergepaths DESTINATION bin)

  # dmri_trk2trk
  add_executable(dmri_trk2trk dmri_trk2trk.cxx spline.cxx vial.cxx TrackIO.cxx)
  target_link_libraries(dmri_trk2trk utils fem_elastic tetgen)
  install(TARGETS dmri_trk2trk DESTINATION bin)

  # dmri_vox2vox
  add_executable(dmri_vox2vox dmri_vox2vox.cxx vial.cxx)
  target_link_libraries(dmri_vox2vox utils fem_elastic tetgen)
  install(TARGETS dmri_vox2vox DESTINATION bin)

  # dmri_train
  add_executable(dmri_train dmri_train.cxx blood.cxx spline.cxx vial.cxx TrackIO.cxx)
  target_link_libraries(dmri_train utils fem_elastic tetgen)
  install(TARGETS dmri_train DESTINATION bin)

  # dmri_paths
  add_executable(dmri_paths dmri_paths.cxx coffin.cxx bite.cxx spline.cxx vial.cxx TrackIO.cxx)
  target_link_libraries(dmri_paths utils fem_elastic tetgen)
  install(TARGETS dmri_paths DESTINATION bin)

  # dmri_pathstats
  add_executable(dmri_pathstats dmri_pathstats.cxx spline.cxx blood.cxx vial.cxx TrackIO.cxx)
  target_link_libraries(dmri_pathstats utils fem_elastic tetgen)
  install(TARGETS dmri_pathstats DESTINATION bin)

  # dmri_motion
  add_executable(dmri_motion dmri_motion.cxx vial.cxx)
  target_link_libraries(dmri_motion utils fem_elastic tetgen)
  install(TARGETS dmri_motion DESTINATION bin)

  # dmri_group
  add_executable(dmri_group dmri_group.cxx vial.cxx)
  target_link_libraries(dmri_group utils fem_elastic tetgen)
  install(TARGETS dmri_group DESTINATION bin)

endif()
