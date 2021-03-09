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

    {
        typedef void FORVAR;    // make sure not shared
        if (sizeof(FORVAR*)) {
	  ;   // use to avoid a warning msg
	}
        
        // Decide the distribution of the loops
        // and start the sum reductions, if any
        //
        ROMP_Distributor ROMP_distributor;

        ROMP_Distributor_begin(&ROMP_distributor,
            ROMP_LO, ROMP_HI
#ifdef ROMP_SUMREDUCTION0
            , &ROMP_SUMREDUCTION0
#else
            , NULL
#endif

#ifdef ROMP_SUMREDUCTION1
            , &ROMP_SUMREDUCTION1
#else
            , NULL
#endif

#ifdef ROMP_SUMREDUCTION2
            , &ROMP_SUMREDUCTION2
#else
            , NULL
#endif
            );

    {
        // Poison the reduced identifiers
#ifdef ROMP_SUMREDUCTION0
      typedef void ROMP_SUMREDUCTION0;
      if (sizeof(ROMP_SUMREDUCTION0*)) {
	;  // use to avoid a warning msg
      }
#endif
#ifdef ROMP_SUMREDUCTION1
      typedef void ROMP_SUMREDUCTION1;
      if (sizeof(ROMP_SUMREDUCTION1*)) {
	;  // use to avoid a warning msg
      }
#endif
#ifdef ROMP_SUMREDUCTION2
      typedef void ROMP_SUMREDUCTION2;
      if (sizeof(ROMP_SUMREDUCTION2*)) {
	;  // use to avoid a warning msg
      }
#endif
                
        // Parallel iteration over the partial sums
        //       
        int ROMP_index;

        ROMP_PF_begin
#ifdef ROMP_SUPPORT_ENABLED
        ROMP_pf_static.line = romp_for_line;
#endif

// When this is not in a macro, but is directly here instead, the compilers use this file 
// as the source correlation of the "omp parallel for" loop body, which hinders using performance analysis tools
//
#ifndef ROMP_for_begin
#define ROMP_for_begin                                                                      \
  	for (ROMP_index = 0; ROMP_index < ROMP_distributor.partialSize; ROMP_index++) {     \
            ROMP_PF_begin                                                                   \
                                                                                            \
            /* Serial iteration reproducing each partial sum    */                          \
            /*                                                  */                          \
            int const ROMP_lo = ROMP_distributor.partials[ROMP_index].lo;                   \
            int const ROMP_hi = ROMP_distributor.partials[ROMP_index].hi;                   \
            int ROMP_VARIABLE;                                                              \
            for (ROMP_VARIABLE = ROMP_lo; ROMP_VARIABLE < ROMP_hi; ROMP_VARIABLE++) {       \
                                                                                            \
	        ROMP_PFLB_begin                                                             \
    // end of macro
#endif

#ifndef ROMP_FOR_LEVEL
#define ROMP_FOR_LEVEL assume_reproducible;
#endif

#ifdef HAVE_OPENMP
#pragma omp parallel for if_ROMPLEVEL(ROMP_FOR_LEVEL)
#endif
// ROMP_for_begin should be the next line after the include
