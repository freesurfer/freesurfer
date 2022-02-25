project(julia)

install(FILES
  julia/FreeSurfer/Manifest.toml
  julia/FreeSurfer/Project.toml
  DESTINATION julia/FreeSurfer
)

install(FILES
  julia/FreeSurfer/src/FreeSurfer.jl
  julia/FreeSurfer/src/dti.jl
  julia/FreeSurfer/src/mri.jl
  julia/FreeSurfer/src/trk.jl
  julia/FreeSurfer/src/view.jl
  DESTINATION julia/FreeSurfer/src
)

