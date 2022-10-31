project(mri_synthmorph)

# The model file archive is in the git-annex. If the annex files
# are not checked out, then the model archive will not be here
# and no model files will be expanded or installed.
set(MODELARCH ${CMAKE_CURRENT_SOURCE_DIR}/mu40-crazy.tgz)
set(MODELDIR ${CMAKE_CURRENT_SOURCE_DIR}/mu40-crazy)
if(EXISTS ${MODELARCH})
   if(NOT EXISTS ${MODELDIR})
      execute_process(COMMAND bash -c "${CMAKE_CURRENT_SOURCE_DIR}/expand_model_files.sh ${CMAKE_CURRENT_SOURCE_DIR} mu40-crazy" RESULT_VARIABLE retcode)
   endif()
endif()

# Always install the script which can get model files form local paths
install_pyscript(mri_synthmorph)

# Only install model files if git-annex file mu40-crazy.tgz has been found and processed above
if(EXISTS ${MODELDIR})
   install(FILES  mu40-crazy/seg_feat_neck_noise_3_3d_256f/10000.h5 DESTINATION models/mu40-crazy/seg_feat_neck_noise_3_3d_256f)
   install(FILES  mu40-crazy/seg_feat_neck_rigid_3d_256f/10000.h5   DESTINATION models/mu40-crazy/seg_feat_neck_rigid_3d_256f)
   install(FILES  mu40-crazy/weights_def_half.h5 DESTINATION models/mu40-crazy)
endif()

