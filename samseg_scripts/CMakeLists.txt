project(samseg_repo_scripts)

if (FSPYTHON_INSTALL_TREE)
   set(pydir "${FSPYTHON_INSTALL_PREFIX}/python")
else()
   set(pydir "${CMAKE_INSTALL_PREFIX}/python")
endif()

# These local copies of samseg repo scripts replace the install of any gems versions (if one exists
# by the same name under ./freessurfer/samseg) when -DINTEGRATE_SAMSEG=ON.  This allows for:
#
# (1) the scripts to be under source control in (and reproduced from) the freesurfer github repo 
#     as part of the standard git checkout and configure/build/install process.  Builds expects to
#     directly reproduce the source tree by by checking out the source files via a git tag.
# (2) use of the existing install functions, e.g., install_pyscript(), install_pyscript_fspython_tree(), etc.
#
# (1) and (2) apply to all scripts installed under python/scripts. The python install functions install
# scripts with the correct privileges there and generate and install a wrapper script under a ./bin subdir 
# using paths for either the single install tree or split install tree (when -DFSPYTHON_INSTALL_TREE=ON).
# However, no fspython package installs files under under ./python/scripts

if(INTEGRATE_SAMSEG)
   if (FSPYTHON_INSTALL_TREE)
      install_pyscript_fspython_tree(computeTissueConcentrations gems_compute_atlas_probs merge_add_mesh_alphas prepareAtlasDirectory run_samseg run_samseg_long sbtiv segment_subregions)
   else()
      install_pyscript(computeTissueConcentrations gems_compute_atlas_probs merge_add_mesh_alphas prepareAtlasDirectory run_samseg run_samseg_long sbtiv segment_subregions)
   endif()
endif()

