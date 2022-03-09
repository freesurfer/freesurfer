#=
  Original Author: Anastasia Yendiki

  Copyright © 2022 The General Hospital Corporation (Boston, MA) "MGH"
 
  Terms and conditions for use, reproduction, distribution and contribution
  are found in the 'FreeSurfer Software License Agreement' contained
  in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 
  https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 
  Reporting: freesurfer@nmr.mgh.harvard.edu
=#

module FreeSurfer

include("mri.jl")
include("trk.jl")
include("view.jl")
include("dti.jl")

function __init__()

  println("FREESURFER_HOME: " * (haskey(ENV, "FREESURFER_HOME") ?
                                 ENV["FREESURFER_HOME"] : "not defined"))
end

end # module
