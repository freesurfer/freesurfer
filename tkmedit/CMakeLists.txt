project(tkmedit)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(tkmedit tkmedit.c)
target_link_libraries(tkmedit utils)

install(TARGETS tkmedit DESTINATION bin)
