project(fsdeepnet_repo_scripts)

if (FSPYTHON_INSTALL_TREE)
   set(pydir "${FSPYTHON_INSTALL_PREFIX}/python")
else()
   set(pydir "${CMAKE_INSTALL_PREFIX}/python")
endif()

if (FSPYTHON_INSTALL_TREE)
   install_pyscript_fspython_tree(fsdeepnet_create_data_list fsdeepnet_predict fsdeepnet_test_dataset fsdeepnet_train fsdeepnet_checkpoint fsdeepnet_evaluate fsdeepnet_test_dataloader fsdeepnet_tf2torch)
else()
   install_pyscript(fsdeepnet_create_data_list fsdeepnet_predict fsdeepnet_test_dataset fsdeepnet_train fsdeepnet_checkpoint fsdeepnet_evaluate fsdeepnet_test_dataloader fsdeepnet_tf2torch)
endif()

