project(dmri_tensoreig)

include_directories(${FS_INCLUDE_DIRS})

add_executable(dmri_tensoreig dmri_tensoreig.c)
target_link_libraries(dmri_tensoreig utils)

install(TARGETS dmri_tensoreig DESTINATION bin)
