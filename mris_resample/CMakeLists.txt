project(mris_resample)

if(ANN_FOUND AND Boost_FOUND)

  include_directories(${FS_INCLUDE_DIRS} ${ANN_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})

  add_executable(mris_resample mris_resample.cpp)
  target_link_libraries(mris_resample utils ${ANN_LIBRARIES} ${Boost_PROGRAM_OPTIONS_LIBRARY})

  install(TARGETS mris_resample DESTINATION bin)

endif()
