project(fsdeepnet_models)

set(FSDEEPNET_MODEL_SERVER "https://ftp.nmr.mgh.harvard.edu/pub/dist/lcnpublic/dist/fsdeepnet-models")
set(FSDEEPNET_PRETRAINED
    "entowm.fsm31.t1.nstd00-30.nstd21-108.pth"
    "mca-dura.both-lh.nstd21.fhs.pth"
    "sclimbic.fsm+ad.t1.nstd00-50.nstd32-50.pth"
    "vsinus.no-sp.m.all.nstd10-070.pth"
    "claustrum_seg_20250616.pth"
    "synthstrip.1.pth"
    "synthstrip.nocsf.1.pth"
    "synthseg_2.0.pth"
    "synthseg_parc_2.0.pth"
    "pglands_seg.pth"
    )

function(download_fsdeepnet_pretrained)
  foreach(PRETRAINED ${ARGN})
    install(CODE "
      message(STATUS \"Downloading fsdeepnet model file: ${CMAKE_INSTALL_PREFIX}/fsdeepnet_models/${PRETRAINED}\")
      file(DOWNLOAD
           \"${FSDEEPNET_MODEL_SERVER}/${PRETRAINED}\"
           \"${CMAKE_INSTALL_PREFIX}/fsdeepnet_models/${PRETRAINED}\"
           STATUS download_status)

      # Check if the download was successful
      list(GET download_status 0 status_code)
      if(NOT status_code EQUAL 0)
         message(FATAL_ERROR \"Download ${CMAKE_INSTALL_PREFIX}/fsdeepnet_models/${PRETRAINED} failed\")
      endif()"
    )
  endforeach()
endfunction()

download_fsdeepnet_pretrained(${FSDEEPNET_PRETRAINED})
