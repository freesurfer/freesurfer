project(ThalamicNucleiDTI)

if(APPLE)
    # TODO Cereate the MacOS executable, add the call to set MacOS BIN_PROGRAMS
    message(STATUS "ThalamicNucleiDTI is not currently available for MacOS, skipping install")

else()
    set(BIN_PROGRAMS
        linux_x86_64/SegmentThalamicNuclei_DTI          
        linux_x86_64/run_SegmentThalamicNuclei_DTI.sh   
    )

    install_symlinks(
        ${BIN_PROGRAMS}
        segmentThalamicNuclei_DTI.sh
        TYPE programs
        DESTINATION bin
    )
endif()

file(GLOB ATLAS_FILES
    atlases/"*")

install_symlinks(
    ${ATLAS_FILES}
    TYPE files
    DESTINATION  average/ThalamicNuclei/atlas_DTI
)
