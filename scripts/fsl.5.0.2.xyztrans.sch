clear UT

## This is an FSL file. It is copied of xyztrans.sch from FSL version 
## 5.0.2. It is used in fslregister

## X TRANSLATION ##

# 8mm scale
setscale 8
setoption smoothing 8
setoption boundguess 8
setoption paramsubset 1  0 0 0 1 0 0 0 0 0 0 0 0
clear U
clear UA
setrow UA 1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1
optimise 12 UA:1   0.0   0.0   0.0    0.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0    8.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   -8.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   16.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0  -16.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   24.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0  -24.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   32.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0  -32.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   40.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0  -40.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   48.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0  -48.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   56.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0  -56.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   64.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0  -64.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   72.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0  -72.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   80.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0  -80.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   88.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0  -88.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   96.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0  -96.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0  104.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0 -104.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0  112.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0 -112.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0  120.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0 -120.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0  128.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0 -128.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0  136.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0 -136.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0  144.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0 -144.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0  152.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0 -152.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0  160.0   0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0 -160.0   0.0   0.0   0.0   abs 4 
clear UA
copy U UA

# 4mm scale
setscale 4
setoption smoothing 4
setoption boundguess 4
setoption paramsubset 1  0 0 0 1 0 0 0 0 0 0 0 0
clear UB
clear UL
clear UM
# remeasure costs at this scale
clear U
measurecost 12 UA 0 0 0 0 0 0 abs
sort U
copy U UL
# optimise best 3 candidates
clear U
optimise 12 UL:1-3  0.0   0.0   0.0   0.0   0.0   0.0   0.0  abs 4
# also try the identity transform as a starting point at this resolution
clear UQ
setrow UQ  1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1
optimise 12 UQ  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
clear UB
copy U UB

# 2mm scale
setscale 2
setoption smoothing 2
setoption boundguess 2
setoption paramsubset 1  0 0 0 1 0 0 0 0 0 0 0 0
clear U
clear UC
clear UF
# remeasure costs at this scale
measurecost 12 UB 0 0 0 0 0 0 abs
sort U
copy U UC
clear U
optimise 12  UC:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  abs 4
copy U UF

# 1mm scale
setscale 1
setoption smoothing 1
setoption boundguess 1
setoption paramsubset 1  0 0 0 1 0 0 0 0 0 0 0 0
clear U
# also try the identity transform as a starting point at this resolution
setrow UF  1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1
optimise 12 UF:1-2  0.0   0.0   0.0   0.0   0.0   0.0   0.0  abs 1
sort U
copy U:1 UT


## Y TRANSLATION ##

# 8mm scale
setscale 8
setoption smoothing 8
setoption boundguess 8
setoption paramsubset 1  0 0 0 0 1 0 0 0 0 0 0 0
clear U
clear UA
setrow UA 1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1
optimise 12 UA:1   0.0   0.0   0.0   0.0    0.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0    8.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0   -8.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0   16.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0  -16.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0   24.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0  -24.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0   32.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0  -32.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0   40.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0  -40.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0   48.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0  -48.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0   56.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0  -56.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0   64.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0  -64.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0   72.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0  -72.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0   80.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0  -80.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0   88.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0  -88.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0   96.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0  -96.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0  104.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0 -104.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0  112.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0 -112.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0  120.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0 -120.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0  128.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0 -128.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0  136.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0 -136.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0  144.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0 -144.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0  152.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0 -152.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0  160.0   0.0   0.0   abs 4 
optimise 12 UA:1   0.0   0.0   0.0   0.0 -160.0   0.0   0.0   abs 4 
clear UA
copy U UA

# 4mm scale
setscale 4
setoption smoothing 4
setoption boundguess 4
setoption paramsubset 1  0 0 0 0 1 0 0 0 0 0 0 0
clear UB
clear UL
clear UM
# remeasure costs at this scale
clear U
measurecost 12 UA 0 0 0 0 0 0 abs
sort U
copy U UL
# optimise best 3 candidates
clear U
optimise 12 UL:1-3  0.0   0.0   0.0   0.0   0.0   0.0   0.0  abs 4
# also try the identity transform as a starting point at this resolution
clear UQ
setrow UQ  1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1
optimise 7 UQ  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
clear UB
copy U UB

# 2mm scale
setscale 2
setoption smoothing 2
setoption boundguess 2
setoption paramsubset 1  0 0 0 0 1 0 0 0 0 0 0 0
clear U
clear UC
clear UF
# remeasure costs at this scale
measurecost 12 UB 0 0 0 0 0 0 abs
sort U
copy U UC
clear U
optimise 12  UC:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  abs 4
copy U UF

# 1mm scale
setscale 1
setoption smoothing 1
setoption boundguess 1
setoption paramsubset 1  0 0 0 0 1 0 0 0 0 0 0 0
clear U
# also try the identity transform as a starting point at this resolution
setrow UF  1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1
optimise 12 UF:1-2  0.0   0.0   0.0   0.0   0.0   0.0   0.0  abs 1
sort U
copy U:1 UT



## Z TRANSLATION ##

# 8mm scale
setscale 8
setoption smoothing 8
setoption boundguess 8
setoption paramsubset 1  0 0 0 0 0 1 0 0 0 0 0 0
clear U
clear UA
setrow UA 1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0    0.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0    8.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0   -8.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0   16.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0  -16.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0   24.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0  -24.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0   32.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0  -32.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0   40.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0  -40.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0   48.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0  -48.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0   56.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0  -56.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0   64.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0  -64.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0   72.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0  -72.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0   80.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0  -80.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0   88.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0  -88.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0   96.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0  -96.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0  104.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0 -104.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0  112.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0 -112.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0  120.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0 -120.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0  128.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0 -128.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0  136.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0 -136.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0  144.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0 -144.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0  152.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0 -152.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0  160.0   0.0   abs 4 
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0 -160.0   0.0   abs 4 
clear UA
copy U UA

# 4mm scale
setscale 4
setoption smoothing 4
setoption boundguess 4
setoption paramsubset 1  0 0 0 0 0 1 0 0 0 0 0 0
clear UB
clear UL
clear UM
# remeasure costs at this scale
clear U
measurecost 12 UA 0 0 0 0 0 0 abs
sort U
copy U UL
# optimise best 3 candidates
clear U
optimise 12 UL:1-3  0.0   0.0   0.0   0.0   0.0   0.0   0.0  abs 4
# also try the identity transform as a starting point at this resolution
clear UQ
setrow UQ  1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1
optimise 7 UQ  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
clear UB
copy U UB

# 2mm scale
setscale 2
setoption smoothing 2
setoption boundguess 2
setoption paramsubset 1  0 0 0 0 0 1 0 0 0 0 0 0
clear U
clear UC
clear UF
# remeasure costs at this scale
measurecost 12 UB 0 0 0 0 0 0 abs
sort U
copy U UC
clear U
optimise 12  UC:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  abs 4
copy U UF

# 1mm scale
setscale 1
setoption smoothing 1
setoption boundguess 1
setoption paramsubset 1  0 0 0 0 0 1 0 0 0 0 0 0
clear U
# also try the identity transform as a starting point at this resolution
setrow UF  1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1
optimise 12 UF:1-2  0.0   0.0   0.0   0.0   0.0   0.0   0.0  abs 1
sort U
copy U:1 UT

## sort the 3 results to pick the best
clear U
copy UT U
sort U

# now do a general 3 DOF translation to refine this
clear UA
copy U UA

# 8mm scale
setscale 8
setoption smoothing 8
setoption paramsubset 3  0 0 0 1 0 0 0 0 0 0 0 0  0 0 0 0 1 0 0 0 0 0 0 0  0 0 0 0 0 1 0 0 0 0 0 0
clear U
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4 

# 4mm scale
setscale 4
setoption smoothing 4
setoption paramsubset 3  0 0 0 1 0 0 0 0 0 0 0 0  0 0 0 0 1 0 0 0 0 0 0 0  0 0 0 0 0 1 0 0 0 0 0 0
clear UB
clear UL
clear UM
# remeasure costs at this scale
clear U
measurecost 12 UA 0 0 0 0 0 0 rel
sort U
copy U UL
# optimise best 3 candidates
clear U
optimise 12 UL:1-3  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
# also try the identity transform as a starting point at this resolution
clear UQ
setrow UQ  1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1
optimise 7 UQ  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
copy U UB

# 2mm scale
setscale 2
setoption smoothing 2
setoption paramsubset 3  0 0 0 1 0 0 0 0 0 0 0 0  0 0 0 0 1 0 0 0 0 0 0 0  0 0 0 0 0 1 0 0 0 0 0 0
clear U
clear UC
clear UD
clear UE
clear UF
# remeasure costs at this scale
measurecost 12 UB 0 0 0 0 0 0 rel
sort U
copy U UC
clear U
optimise 12  UC:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
sort U
copy U UF

# 1mm scale
setscale 1
setoption smoothing 1
setoption boundguess 1
setoption paramsubset 3  0 0 0 1 0 0 0 0 0 0 0 0  0 0 0 0 1 0 0 0 0 0 0 0  0 0 0 0 0 1 0 0 0 0 0 0
clear U
# also try the identity transform as a starting point at this resolution
setrow UF  1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1
optimise 12 UF:1-2  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 1
sort U

