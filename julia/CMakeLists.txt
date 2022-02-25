project(julia)

install(FILES
  FreeSurfer/Manifest.toml
  FreeSurfer/Project.toml
  DESTINATION julia/FreeSurfer
)

install(FILES
  FreeSurfer/src/FreeSurfer.jl
  FreeSurfer/src/dti.jl
  FreeSurfer/src/mri.jl
  FreeSurfer/src/trk.jl
  FreeSurfer/src/view.jl
  DESTINATION julia/FreeSurfer/src
)

