#
# This is a schedule file given to us (FreeSurfer, DNG) on 8/11/08 to
# fix the "intialization" problem.  Here's MJ's msg:
#
#   It turns out that it is largely to do with the 8mm phase not being
#   particularly good for searching the cost function.
#   I have fixed the problem by explicitly trying an identity starting
#   matrix within the schedule (this identity is taken in conjunction
#   with the -init) at both 4mm and 1mm phases.  This seems to fix the
#   problem completely. I have now made this change in our code base
#   and it should get shipped with the next fsl release
#

# 8mm scale
setscale 8
setoption smoothing 8
clear S
clear P
search
# 4mm scale
setscale 4
setoption smoothing 4
clear U
clear UA 
clear UB
clear US
clear UP
# remeasure costs at this scale
measurecost 7 S 0 0 0 0 0 0 rel
copy U US
clear U
measurecost 7 P 0 0 0 0 0 0 rel
copy U UP
dualsort US UP
# optimise best 3 candidates (pre and post 8mm optimisations)
clear U
optimise 7 US:1-3  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
optimise 7 UP:1-3  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
# also try the identity transform as a starting point at this resolution
clear UQ
setrow UQ  1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1
optimise 7 UQ  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
sort U
copy U UA
# select best 4 optimised solutions and try perturbations of these
clear U
copy UA:1-4 U
optimise 7 UA:1-4  1.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
optimise 7 UA:1-4 -1.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
optimise 7 UA:1-4  0.0   1.0   0.0   0.0   0.0   0.0   0.0  rel 4
optimise 7 UA:1-4  0.0  -1.0   0.0   0.0   0.0   0.0   0.0  rel 4
optimise 7 UA:1-4  0.0   0.0   1.0   0.0   0.0   0.0   0.0  rel 4
optimise 7 UA:1-4  0.0   0.0  -1.0   0.0   0.0   0.0   0.0  rel 4
optimise 7 UA:1-4  0.0   0.0   0.0   0.0   0.0   0.0   0.1  abs 4
optimise 7 UA:1-4  0.0   0.0   0.0   0.0   0.0   0.0  -0.1  abs 4
optimise 7 UA:1-4  0.0   0.0   0.0   0.0   0.0   0.0   0.2  abs 4
optimise 7 UA:1-4  0.0   0.0   0.0   0.0   0.0   0.0  -0.2  abs 4
sort U
copy U UB
# 2mm scale
setscale 2
setoption smoothing 2
clear U
clear UC
clear UD
clear UE
clear UF
# remeasure costs at this scale
measurecost 7 UB 0 0 0 0 0 0 rel
sort U
copy U UC
clear U
optimise 7  UC:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
copy U UD
setoption boundguess 1
if MAXDOF > 7
 clear U
if MAXDOF > 7
 optimise 9  UD:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 1
copy U UE
if MAXDOF > 9
 clear U
if MAXDOF > 9
 optimise 12 UE:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 2
sort U
copy U UF
# 1mm scale
setscale 1
setoption smoothing 1
setoption boundguess 1
clear U
# also try the identity transform as a starting point at this resolution
setrow UF  1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1
optimise 12 UF:1-2  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 1
sort U
