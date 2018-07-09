# install all licenses to the doc directory
file(GLOB LICENSES "license*txt")
install(FILES ${LICENSES} DESTINATION docs)
