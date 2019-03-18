project(mkxsubjreg)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mkxsubjreg mkxsubjreg.cpp)
target_link_libraries(mkxsubjreg utils)

install(TARGETS mkxsubjreg DESTINATION bin)
