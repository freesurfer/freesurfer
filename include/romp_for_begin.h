/**
 * @file  romp_support.h
 * @brief prototypes and structures for getting reprodiucible results from and for timing omp loops.
 *
 */
/*
 * Original Author: Bevin Brett
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2017/12 $
 *    $Revision: 1.0 $
 *
 * Copyright © 2012 The General Hospital Corporation (Boston, MA) "MGH"
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
#endif
#ifdef ROMP_SUMREDUCTION1
        typedef void ROMP_SUMREDUCTION1;
#endif
#ifdef ROMP_SUMREDUCTION2
        typedef void ROMP_SUMREDUCTION2;
#endif
                
        // Parallel iteration over the partial sums
        //       
        int ROMP_index;

        ROMP_PF_begin
#ifdef HAVE_OPENMP
#ifdef ROMP_FOR_LEVEL 
    #pragma omp parallel for if_ROMPLEVEL(ROMP_FOR_LEVEL)
#else
    shown_reproducible
    #pragma omp parallel for if_ROMP(assume_reproducible)
#endif
#endif
  	for (ROMP_index = 0; ROMP_index < ROMP_distributor.partialSize; ROMP_index++) {
            ROMP_PFLB_begin
            
            // Serial iteration reproducing each partial sum
            //
            int const ROMP_lo = ROMP_distributor.partials[ROMP_index].lo; 
            int const ROMP_hi = ROMP_distributor.partials[ROMP_index].hi; 
            int ROMP_VARIABLE;
            for (ROMP_VARIABLE = ROMP_lo; ROMP_VARIABLE < ROMP_hi; ROMP_VARIABLE++) {

	        ROMP_PFLB_begin
