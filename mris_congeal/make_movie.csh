#!/bin/tcsh -ef

source subjects.csh
setenv SUBJECTS_DIR $ZILLES

set avg=fsaverage

setenv HEMI lh

setenv LABEL MT
foreach sno (1)
    foreach iter (0000 0001 0002 )
        setenv ITER ${iter}
        setenv CNAME ./MT.congealed.template.sno${sno}.iter${iter}
        setenv SNAME sphere
        setenv SUBJECT fsaverage
        setenv SUFFIX MT.congealed.target.sno${sno}
        tksurfer $SUBJECT $HEMI $SNAME -tcl render_surface.tcl
        setenv CNAME sulc
        setenv SNAME MT.congealed.sno0.iter0.${iter}
        setenv SUFFIX MT.congealed
        foreach s ($subjects)
            setenv SUBJECT $s
#            tksurfer $SUBJECT $HEMI $SNAME -tcl render_surface.tcl
        end
    end
end
