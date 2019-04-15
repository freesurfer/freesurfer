project(mris_resample)

if(ANN_FOUND)
  include_directories(${FS_INCLUDE_DIRS} ${ANN_INCLUDE_DIR})

  add_executable(mris_resample mris_resample.cpp)
  add_help(mris_resample mris_resample.help.xml)
  target_link_libraries(mris_resample utils ${ANN_LIBRARIES})

  install(TARGETS mris_resample DESTINATION bin)
endif()
