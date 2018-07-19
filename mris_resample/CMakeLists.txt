project(mris_resample)

if(Boost_FOUND)

  include_directories(${CMAKE_SOURCE_DIR}/include ${ANN_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})

  add_executable(mris_resample mris_resample.cpp)
  target_link_libraries(mris_resample utils ${ANN_LIBRARIES} ${Boost_PROGRAM_OPTIONS_LIBRARY})

  install(TARGETS mris_resample DESTINATION bin)

endif()
