#=
  Original Author: Anastasia Yendiki

  Copyright © 2022 The General Hospital Corporation (Boston, MA) "MGH"

  Terms and conditions for use, reproduction, distribution and contribution
  are found in the 'FreeSurfer Software License Agreement' contained
  in the file 'LICENSE' found in the FreeSurfer distribution, and here:

  https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense

  Reporting: freesurfer@nmr.mgh.harvard.edu
=#

using ColorTypes, DelimitedFiles, ImageInTerminal

export LUT, color_lut, show, view

"Container for segmentation and tract look-up tables"
mutable struct LUT
  id::Vector{Int}
  name::Vector{String}
  rgb::Vector{RGB}
end


"""
    LUT()

Return an empty `LUT` structure
"""
LUT() = LUT(
  Vector{Int}(undef, 0),
  Vector{String}(undef, 0),
  Vector{RGB}(undef, 0)
)


"""
    LUT(infile::String)

Read a look-up table from `infile` and return a `LUT` structure

The input file is assumed to have the format of FreeSurferColorLUT.txt
"""
function LUT(infile::String)

  lut = LUT()

  if ~isfile(infile)
    error(infile * "is not a regular file")
  end

  tab = readdlm(infile; comments=true, comment_char='#')

  lut.id   = tab[:,1]
  lut.name = tab[:,2]
  lut.rgb  = RGB.(tab[:,3]/255, tab[:,4]/255, tab[:,5]/255)

  return lut
end


"The FreeSurfer color look-up table"
const global color_lut = haskey(ENV, "FREESURFER_HOME") ?
                         LUT(ENV["FREESURFER_HOME"]*"/FreeSurferColorLUT.txt") :
                         LUT()


"""
    vol_to_rgb(vol::Array)

Convert an image array to an RGB/Gray array for display.

Determine how the image should be displayed:
- If all the values are IDs in the FreeSurfer color look-up table,
  the image is assumed to be is a segmentation map and is converted
  to RGB based on the FreeSurfer color look-up table.
- If the image has size 3 in any dimension, and the sum of squares
  of the values in that dimension is approx. 1 everywhere, the image
  is assumed to be a vector map and is converted to RGB based on
  vector orientation.
- Otherwise, the image is assumed to be a generic intensity map and
  is converted to grayscale.

Return an array of RGB/Gray values the same size as the input vol.
"""
function vol_to_rgb(vol::Array)

  if isempty(color_lut.id)
    error("FreeSurfer color look-up table is undefined")
  end

  if ~any(isnothing.(indexin(unique(vol), color_lut.id)))
    # Assume the input is a segmentation map, get RGB of labels from LUT
    return color_lut.rgb[indexin(vol, color_lut.id)]
  end

  dim3 = findall(size(vol) .== 3)

  for idim in dim3
    if all(isapprox.(sum(vol.^2; dims=idim), 1) .|| all(vol.==0, dims=idim))
      # Assume the input is a vector map, get RGB based on vector orientation
      return RGB.(abs.(selectdim(vol,idim,1)),
                  abs.(selectdim(vol,idim,2)),
                  abs.(selectdim(vol,idim,3)))
    end
  end

  # Otherwise, assume the input is a generic intensity map
  return Gray.(vol / maximum(vol))
end


"""
    show(mri::MRI, mrimod::Union{MRI, Nothing}=nothing)

Show an `MRI` structure (an image slice and a summary of header info) in the
terminal window

# Arguments:
- `mri::MRI`: the main image to display
- `mrimod:MRI`: an optional image to modulate the main image by (e.g., an FA
  map to modulate a vector map)
"""
function show(mri::MRI, mrimod::Union{MRI, Nothing}=nothing)

  # Find non-empty slices in z dimension
  iz = any(mri.vol .!= 0; dims=([1:2; 4:ndims(mri.vol)]))
  iz = reshape(iz, length(iz))
  iz = findall(iz)

  # Keep only the middle non-empty slice in z dimension
  iz = iz[Int(round(end/2))]

  # Find non-empty slices in x, y dimensions
  ix = any(mri.vol[:,:,iz,:] .!= 0; dims=(2:ndims(mri.vol)))
  ix = reshape(ix, length(ix))
  ix = findall(ix)

  iy = any(mri.vol[:,:,iz,:] .!= 0; dims=([1; 3:ndims(mri.vol)]))
  iy = reshape(iy, length(iy))
  iy = findall(iy)

  # Keep only the area containing non-empty slices in x, y dimensions
  ix = ix[1]:ix[end]
  iy = iy[1]:iy[end]

  # Subsample to fit display in x, y dimensions
  nsub = Int(ceil(length(iy) ./ displaysize()[2]))

  ix = ix[1:nsub:end]
  iy = iy[1:nsub:end]

  # Convert image to RGB/Gray array
  rgb = vol_to_rgb(mri.vol[ix, iy, iz, :])

  # Keep first frame only
  rgb = rgb[:, :, 1]

  # Optional intensity modulation
  if ~isnothing(mrimod)
    if size(mrimod.vol)[1:3] != size(mri.vol)[1:3]
      error("Dimension mismatch between main image " * 
            string(size(mri.vol)[1:3]) *
            " and modulation image " *
            string(size(mrimod.vol)[1:3]))
    end

    rgbmod = mrimod.vol[ix, iy, iz, 1] / maximum(mrimod.vol)
  end

  if hasfield(eltype(rgb), :r)		# RGB
    # Add α channel to make zero voxels transparent
    rgb = RGBA.(getfield.(rgb, :r), getfield.(rgb, :g), getfield.(rgb, :b),
                Float64.(rgb .!= RGB(0,0,0)))

    # Modulate intensity with (optional) second image
    if ~isnothing(mrimod)
      rgb = RGBA.(getfield.(rgb, :r) .* rgbmod,
                  getfield.(rgb, :g) .* rgbmod,
                  getfield.(rgb, :b) .* rgbmod,
                  getfield.(rgb, :alpha))
    end
  else					# Gray
    # Add α channel to make zero voxels transparent
    rgb = GrayA.(getfield.(rgb, :val), Float64.(rgb .!= Gray(0)))

    # Modulate intensity with (optional) second image
    if ~isnothing(mrimod)
      rgb = GrayA.(getfield.(rgb, :val) .* rgbmod,
                   getfield.(rgb, :alpha))
    end
  end

  # Display image
  if mri.ispermuted
    imshow(rgb)
  else
    imshow(permutedims(rgb, [2; 1]))
  end

  # Print image info
  println()
  if ~isempty(mri.fspec)
    println("Read from: " * mri.fspec)
  end
  println("Volume dimensions: " * string(collect(size(mri.vol))))
  println("Spatial resolution: " * string(Float64.(mri.volres)))
  if ~isempty(mri.bval)
    println("b-values: " * string(Float64.(unique(mri.bval))))
  end
  println("Intensity range: " * string(Float64.([minimum(mri.vol),
                                                 maximum(mri.vol)])))
end


"""
    view(mri::MRI)

View an `MRI` structure in a slice viewer
"""
function view(mri::MRI)

#=
=#

end

