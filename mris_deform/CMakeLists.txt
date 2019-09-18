project(mris_deform)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_deform mris_deform.cpp)
target_link_libraries(mris_deform utils)
install(TARGETS mris_deform DESTINATION bin)

add_executable(mris_ca_deform mris_ca_deform.cpp)
target_link_libraries(mris_ca_deform utils)
install(TARGETS mris_ca_deform DESTINATION bin)

add_executable(mris_nudge mris_nudge.cpp)
target_link_libraries(mris_nudge utils)
install(TARGETS mris_nudge DESTINATION bin)
