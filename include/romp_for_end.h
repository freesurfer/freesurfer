/**
 * @brief prototypes and structures for getting reprodiucible results from and for timing omp loops.
 *
 */
/*
 * Original Author: Bevin Brett
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


                ROMP_PFLB_end
	    }           // serial loop
        
            ROMP_PF_end
        }               // parallel loop
        ROMP_PF_end
        
    } // poison ROMP_SUMREDUCTION0 etc.
    
#undef ROMP_FOR_PRAGMA
#undef ROMP_VARIABLE
#undef ROMP_LO
#undef ROMP_HI

#ifdef ROMP_FOR_LEVEL
#undef ROMP_FOR_LEVEL
#endif

#ifdef ROMP_SUMREDUCTION0
#undef ROMP_SUMREDUCTION0
#endif

#ifdef ROMP_SUMREDUCTION1
#undef ROMP_SUMREDUCTION1
#endif

#ifdef ROMP_SUMREDUCTION2
#undef ROMP_SUMREDUCTION2
#endif

       ROMP_Distributor_end(&ROMP_distributor);
    }
