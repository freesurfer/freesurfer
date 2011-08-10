usage:

    [program] flow [videoFileName] [frameCount]
    [program] pv [videoFileName] [frameCount]

windows build:

    - use the ParticleVideoSample.sln project with Visual Studio 2008 / 2010
    - you will need to set these environment variables:
        - SBL: path to SimpleBaseLib
        - PVL: path to ParticleVideoLib
        - OPENCV: path to OpenCV-2.2.0 files
        - WINSDK: path to Windows SDK files
	  (e.g., C:\Program Files\Microsoft SDKs\Windows\v6.0A)

linux build:

    - you will need to set these environment variables:
        - SBL: path to SimpleBaseLib
        - PVL: path to ParticleVideoLib
        - OPENCV: path to OpenCV-2.2.0 files
    - first build SBL (run make in SimpleBaseLib/build)
    - then build PVL (run make in ParticleVideoLib/build)
    - then build the sample (run make in ParticleVideoLib/sample)
    - feel free to change/remove the architecture flags

notes:

    - program path should include path.conf that specifies dataPath
    - input files should be in dataPath
    - simpleParticle.conf and varMotion.conf should be in dataPath
    - pv command expects flow to be estimated first
    - a flow frameCount of 2 means a single frame of flow will 
        be computed between the first and second frames
    - the particle clustering (enabled by setting doClustering to true) currently works on in Windows
        (if you know of a fix for Linux, let me know; I haven't investigated it yet)
    - to generate PVL documentation with doxygen, run "doxygen pvl_doxygen.txt" from the doc directory

output:

    - both flow and particle estimation will save diagnostic videos
    - particle set will be saved in simple text format (see below)
    - set saveText in estimateOpticalFlow to true if you want to 
        export flow fields in the simple text format (see below)

particle set file simplified text format:

    - first line is number of particles
    - each additional line is one particle
    - particle line consists of startFrame, endFrame, point count, and sequence of x, y points

motion field file simplified text format:

    - first line gives width and height
    - each additional line gives one row of pixels (bottom to top)
    - each pixel has a u value, v value, and occlusion value (0-255)
