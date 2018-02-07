static int MRIScomputeBorderValues_new(
    MRI_SURFACE * const mris,
    MRI         * const mri_brain,
    MRI         * const mri_smooth,
    double        const inside_hi,
    double        const border_hi,
    double        const border_low,
    double        const outside_low,
    double        const outside_hi,
    double        const sigma,
    float         const max_thickness,
    FILE        * const log_fp,
    int           const which,
    MRI *         const mri_mask,
    double        const thresh,
    int           const flags,
    MRI *         const mri_aseg) 
{
  float const step_size = mri_brain->xsize / 2;
    

  MRI *mri_tmp;
  if (mri_brain->type == MRI_UCHAR) {
    mri_tmp = MRIreplaceValues(mri_brain, NULL, 255, 0);
  } else {
    mri_tmp = MRIcopy(mri_brain, NULL);
  }

  MRISclearMarks(mris); /* for soap bubble smoothing later */


  // Various sums
  //
  float mean_dist = 0, mean_in = 0, mean_out = 0, mean_border = 0;

  // Various counters
  //
  int total_vertices = 0;
  int ngrad_max = 0, ngrad = 0, nmin = 0;
  int nmissing  = 0, nout  = 0, nin  = 0, nfound = 0, nalways_missing = 0, num_changed = 0;
  
  double max_mag_dist = 0.0f;

  // Prepare to map all the surface points to voxels
  //
  MRIS_SurfRAS2VoxelMap* sras2v_map = 
    MRIS_makeRAS2VoxelMap(mri_brain, mris);
  
  // OLD CODE HAS THIS COMMENT, DONT KNOW WHY
  // first compute intensity of local gray/white boundary

  int vno;

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(serial) \
    reduction(+:mean_dist,mean_in,mean_out,mean_border) \
    reduction(+:total_vertices,ngrad_max,ngrad,nmin,nmissing,nout,nin,nfound,nalways_missing,num_changed)
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin
    
    VERTEX * const v = &mris->vertices[vno];
    if (v->ripflag) {
      ROMP_PF_continue;
    }

    if (vno == Gdiag_no) {
      DiagBreak();
    }

    // calculate the unit-length normal to the vertex in voxel space
    // 
    float nx,ny,nz;
    {
      double x,y,z;
      
      double xw, yw, zw;
      x = v->x;
      y = v->y;
      z = v->z;
      MRIS_useRAS2VoxelMap(sras2v_map, mri_brain,x, y, z, &xw, &yw, &zw);
      
      double xw1, yw1, zw1;
      x = v->x + v->nx;
      y = v->y + v->ny;
      z = v->z + v->nz;
      MRIS_useRAS2VoxelMap(sras2v_map, mri_brain,x, y, z, &xw1, &yw1, &zw1);
    
      nx = xw1 - xw;
      ny = yw1 - yw;
      nz = zw1 - zw;

      float dist = sqrt(SQR(nx) + SQR(ny) + SQR(nz));
      if (FZERO(dist)) ROMP_PF_continue;                                            // WAS "dist = 1;" BUT THAT MAKES NO SENSE
      nx /= dist;
      ny /= dist;
      nz /= dist;
    }
    
    /*
      find the distance in the directions parallel and anti-parallel to
      the surface normal in which the gradient is pointing 'inwards'.
      The border will then be constrained to be within that region.
    */
    double inward_dist  =  1.0;
    double outward_dist = -1.0;

    double current_sigma;
    for (current_sigma = sigma; current_sigma <= 10 * sigma; current_sigma *= 2) {
    
      double mag = -1.0;
      float dist;
      for (dist = 0; dist > -max_thickness; dist -= step_size) {

        double dx = v->x - v->origx;
        double dy = v->y - v->origy;
        double dz = v->z - v->origz;
        double orig_dist = fabs(dx * v->nx + dy * v->ny + dz * v->nz);              // THIS FORMULA IS WRONG, SINCE +ve and -ve terms should not cancel
                                                                                    // SHOULD BE sqrt(SQR(...) + SQR(...) + SQR(...))

        if (fabs(dist) + orig_dist > max_thickness) {
          break;
        }

        double xw, yw, zw;
        double const x = v->x + v->nx * dist;
        double const y = v->y + v->ny * dist;
        double const z = v->z + v->nz * dist;
        MRIS_useRAS2VoxelMap(sras2v_map, mri_brain,x, y, z, &xw, &yw, &zw);
        
        MRIsampleVolumeDerivativeScale(mri_tmp, xw, yw, zw, nx, ny, nz, &mag, current_sigma);   // expensive
        if (mag >= 0.0) {
          break;
        }
        
        double val;
        MRIsampleVolume(mri_brain, xw, yw, zw, &val);
        if (val > border_hi) {
          break;
        }
        if (mri_mask) {
          MRIsampleVolume(mri_mask, xw, yw, zw, &val);
          if (val > thresh) {
            break;
          }
        }
      }

      inward_dist = dist + step_size / 2;

      if (DIAG_VERBOSE_ON && mri_brain->xsize < .95 && mag >= 0.0)  // refine inward_dist for hires volumes
      {
        for (dist = inward_dist; dist > -max_thickness; dist -= step_size / 2) {
          double x,y,z;
          
          double xw, yw, zw;

          x = v->x + v->nx * dist;
          y = v->y + v->ny * dist;
          z = v->z + v->nz * dist;
          MRIS_useRAS2VoxelMap(sras2v_map, mri_brain,x, y, z, &xw, &yw, &zw);
          
          double val;
          MRIsampleVolume(mri_brain, xw, yw, zw, &val);

          x = v->x + v->nx * (dist + step_size / 2);
          y = v->y + v->ny * (dist + step_size / 2);
          z = v->z + v->nz * (dist + step_size / 2);
          MRIS_useRAS2VoxelMap(sras2v_map, mri_brain,x, y, z, &xw, &yw, &zw);
          
          double next_val;
          MRIsampleVolume(mri_brain, xw, yw, zw, &next_val);
          
          if (next_val < val)  // found max inwards intensity
          {
            break;
          }
        }
        inward_dist = dist;
      }

      for (dist = 0; dist < max_thickness; dist += step_size) {
        double dx = v->x - v->origx;
        double dy = v->y - v->origy;
        double dz = v->z - v->origz;
        double orig_dist = fabs(dx * v->nx + dy * v->ny + dz * v->nz);
        if (fabs(dist) + orig_dist > max_thickness) {
          break;
        }

        double xw, yw, zw;
        double const x = v->x + v->nx * dist;
        double const y = v->y + v->ny * dist;
        double const z = v->z + v->nz * dist;
        MRIS_useRAS2VoxelMap(sras2v_map, mri_brain,x, y, z, &xw, &yw, &zw);
        MRIsampleVolumeDerivativeScale(mri_tmp, xw, yw, zw, nx, ny, nz, &mag, current_sigma);
        if (mag >= 0.0) {
          break;
        }

        double val;
        MRIsampleVolume(mri_brain, xw, yw, zw, &val);
        if (val < border_low) {
          break;
        }
        if (mri_mask) {
          MRIsampleVolume(mri_mask, xw, yw, zw, &val);
          if (val > thresh) {
            break;
          }
        }
      }

      outward_dist = dist - step_size / 2;
      
      if (!isfinite(outward_dist)) {
        DiagBreak();
      }
      if (inward_dist <= 0 || outward_dist >= 0) {
        break;
      }
    }

    if (inward_dist > 0 && outward_dist < 0) {
      current_sigma = sigma; /* couldn't find anything */
    }

    FILE *fp = NULL;
    if (vno == Gdiag_no) {
      char fname[STRLEN];
      sprintf(fname, "v%d.%2.0f.log", Gdiag_no, sigma * 100);
      fp = fopen(fname, "w");
      fprintf(stdout,
              "v %d: inward dist %2.2f, outward dist %2.2f, sigma %2.1f\n",
              vno,
              inward_dist,
              outward_dist,
              current_sigma);
    }

    v->val2 = current_sigma;
    /*
      search outwards and inwards and find the local gradient maximum
      at a location with a reasonable MR intensity value. This will
      be the location of the edge.
    */

    /* search in the normal direction to find the min value */
    double max_mag_val = -10.0f;
    double max_mag = 0.0f;
    double min_val = 10000.0;
    double min_val_dist = 0.0f;
    int    local_max_found = 0;
    
    float dists[MAX_SAMPLES], mri[MAX_SAMPLES];
    int   numberOfSamples = 0;

    float dist;    
    for (dist = inward_dist; dist <= outward_dist; dist += STEP_SIZE) {

      double val;
      {
        double const x = v->x + v->nx * dist;
        double const y = v->y + v->ny * dist;
        double const z = v->z + v->nz * dist;

        double xw, yw, zw;
        MRIS_useRAS2VoxelMap(sras2v_map, mri_brain,x, y, z, &xw, &yw, &zw);
      
        MRIsampleVolume(mri_brain, xw, yw, zw, &val);
      }
      
      dists[numberOfSamples] = dist;
      mri[numberOfSamples]   = val;
      numberOfSamples++;

      double previous_val;
#if 1
      {
        double const x = v->x + v->nx * (dist - STEP_SIZE);
        double const y = v->y + v->ny * (dist - STEP_SIZE);
        double const z = v->z + v->nz * (dist - STEP_SIZE);
        double xw,yw,zw;
        MRIS_useRAS2VoxelMap(sras2v_map, mri_brain,x, y, z, &xw, &yw, &zw);
        MRIsampleVolume(mri_brain, xw, yw, zw, &previous_val);
      }
#else
      /* find max val within 1 mm in inwards direction */
      {
        previous_val = 0.0;
        float d;
        for (d = 0.25; d <= 1.5; d += 0.25) {
          double xw, yw, zw;

          double const x = v->x + v->nx * (d - 1);
          double const y = v->y + v->ny * (d - 1);
          double const z = v->z + v->nz * (d - 1);
          MRIS_useRAS2VoxelMap(sras2v_map, mri_brain,x, y, z, &xw, &yw, &zw);

          double tmp_val;
          MRIsampleVolume(mri_brain, xw, yw, zw, &tmp_val);

          if (tmp_val > previous_val) {
            previous_val = tmp_val;
          }
        }
      }
#endif

      /* the previous point was inside the surface */
      if (previous_val < inside_hi && previous_val >= border_low) {
      
        double xw, yw, zw;
        double x,y,z;
        double val;
        
        /* see if we are at a local maximum in the gradient magnitude */
        x = v->x + v->nx * dist;
        y = v->y + v->ny * dist;
        z = v->z + v->nz * dist;
        MRIS_useRAS2VoxelMap(sras2v_map, mri_brain,x, y, z, &xw, &yw, &zw);
        MRIsampleVolume(mri_brain, xw, yw, zw, &val);

        x = v->x + v->nx * (dist + STEP_SIZE);
        y = v->y + v->ny * (dist + STEP_SIZE);
        z = v->z + v->nz * (dist + STEP_SIZE);
        MRIS_useRAS2VoxelMap(sras2v_map, mri_brain,x, y, z, &xw, &yw, &zw);
        
        double next_mag;
        MRIsampleVolumeDerivativeScale(mri_tmp, xw, yw, zw, nx, ny, nz, &next_mag, sigma);

        x = v->x + v->nx * (dist - STEP_SIZE);
        y = v->y + v->ny * (dist - STEP_SIZE);
        z = v->z + v->nz * (dist - STEP_SIZE);
        MRIS_useRAS2VoxelMap(sras2v_map, mri_brain,x, y, z, &xw, &yw, &zw);
        
        double previous_mag;
        MRIsampleVolumeDerivativeScale(mri_tmp, xw, yw, zw, nx, ny, nz, &previous_mag, sigma);

        if (val < min_val) {
          min_val = val; /* used if no gradient max is found */
          min_val_dist = dist;
        }

        /* if gradient is big and val is in right range */
        x = v->x + v->nx * dist;
        y = v->y + v->ny * dist;
        z = v->z + v->nz * dist;
        MRIS_useRAS2VoxelMap(sras2v_map, mri_brain,x, y, z, &xw, &yw, &zw);

        double mag;
        MRIsampleVolumeDerivativeScale(mri_tmp, xw, yw, zw, nx, ny, nz, &mag, sigma);
        
        // only for hires volumes - if intensities are increasing don't keep going - in gm
        if ((which == GRAY_WHITE)
        &&  (mri_brain->xsize < .95 || flags & IPFLAG_FIND_FIRST_WM_PEAK) 
        &&  (val > previous_val )
        /* && (next_val > val) */   // NOTE - the old code has this uncommented, but fails to init next_val on many of the paths leading to it!
            ) { 
          break;
        }
 
        if ((mri_aseg != NULL) && (MRIindexNotInVolume(mri_aseg, xw, yw, zw) == 0)) {

          int const label = MRIgetVoxVal(mri_aseg, nint(xw), nint(yw), nint(zw), 0);

          if (vno == Gdiag_no)
            printf("v %d: label distance %2.2f = %s @ (%d %d %d)\n",
                   vno,
                   dist,
                   cma_label_to_name(label),
                   nint(xw),
                   nint(yw),
                   nint(zw));
                   
          if ((mris->hemisphere == LEFT_HEMISPHERE  && IS_RH_CLASS(label)) ||
              (mris->hemisphere == RIGHT_HEMISPHERE && IS_LH_CLASS(label))
             ) {
            if (vno == Gdiag_no)
              printf("v %d: terminating search at distance %2.2f due to presence of contra tissue (%s)\n",
                     vno,
                     dist,
                     cma_label_to_name(label));
            break;
          }
        }
        
        if (which == GRAY_CSF) {
          /*
            sample the next val we would process.
            If it is too low, then we
            have definitely reached the border,
            and the current gradient
            should be considered a local max.

            Don't want to do this for gray/white,
            as the gray/white gradient
            often continues seemlessly into the gray/csf.
          */
          double xw,yw,zw;
          
          double const x = v->x + v->nx * (dist + STEP_SIZE);
          double const y = v->y + v->ny * (dist + STEP_SIZE);
          double const z = v->z + v->nz * (dist + STEP_SIZE);
          MRIS_useRAS2VoxelMap(sras2v_map, mri_brain,x, y, z, &xw, &yw, &zw);
          
          double next_val;
          MRIsampleVolume(mri_brain, xw, yw, zw, &next_val);
          if (next_val < border_low) {
            next_mag = 0;
          }
        }

        if (vno == Gdiag_no) fprintf(fp, "%2.3f  %2.3f  %2.3f  %2.3f  %2.3f\n", dist, val, mag, previous_mag, next_mag);

        /*
          if no local max has been found, or this one
          has a greater magnitude,
          and it is in the right intensity range....
        */
        if (
            /* (!local_max_found || (fabs(mag) > max_mag)) && */
            (fabs(mag) > fabs(previous_mag)) && (fabs(mag) > fabs(next_mag)) && (val <= border_hi) &&
            (val >= border_low)) {
            
          double xw,yw,zw;
          double const x = v->x + v->nx * (dist + 1);
          double const y = v->y + v->ny * (dist + 1);
          double const z = v->z + v->nz * (dist + 1);
          MRIS_useRAS2VoxelMap(sras2v_map, mri_brain,x, y, z, &xw, &yw, &zw);
          
          double next_val;
          MRIsampleVolume(mri_brain, xw, yw, zw, &next_val);
          /*
            if next val is in the right range, and the intensity at
            this local max is less than the one at the previous local
            max, assume it is the correct one.
          */
          if ((next_val >= outside_low) &&
              (next_val <= border_hi  ) &&
              (next_val <= outside_hi ) &&
              (!local_max_found || (max_mag < fabs(mag))))
          {                             // beware, this is non-deterministic! if the max mag has equal fabs(), any could be chosen
            local_max_found = 1;
            max_mag_dist = dist; 
            max_mag      = fabs(mag);
            max_mag_val  = val;         
          }
        }
        else {
          /*
            if no local max found yet, just used largest gradient
            if the intensity is in the right range.
          */
          if ((local_max_found == 0) && (fabs(mag) > max_mag) && (val <= border_hi) && (val >= border_low)) {
            double const x = v->x + v->nx * (dist + 1);
            double const y = v->y + v->ny * (dist + 1);
            double const z = v->z + v->nz * (dist + 1);
            double xw,yw,zw;
            MRIS_useRAS2VoxelMap(sras2v_map, mri_brain,x, y, z, &xw, &yw, &zw);
            
            double next_val;
            MRIsampleVolume(mri_brain, xw, yw, zw, &next_val);
            
            if (next_val >= outside_low && next_val <= border_hi && next_val < outside_hi) {
              max_mag_dist = dist;
              max_mag = fabs(mag);
              max_mag_val = val;
            }
          }
        }
      }
    }

    if (vno == Gdiag_no) {
      fclose(fp);
      fp = NULL;
    }

    // doesn't apply to standard stream - only highres or if user specifies
    if (mri_brain->xsize < .95 || flags & IPFLAG_FIND_FIRST_WM_PEAK) {
#ifdef WSIZE
#undef WSIZE
#endif
#define WSIZE 7
      int const whalf = WSIZE;

      if (vno == Gdiag_no) DiagBreak();
      {
        int n;
        for (n = 0; n < v->vnum; n++)
          if (v->v[n] == Gdiag_no) {
            DiagBreak();
          }
      }
      
      // find max in range, and also compute derivative and put it in dm array
      float max_mri = 0;
      float dm[MAX_SAMPLES];
      {
        int i;
        for (i = 0; i < numberOfSamples; i++) {
          if (mri[i] > max_mri) max_mri = mri[i];
          if (i < numberOfSamples - 1 && i > 0) {
            dm[i] = mri[i + 1] - mri[i - 1];
          } else {
            dm[i] = 0;
          }
        }
      }
      
      // compute second derivative
      float dm2[MAX_SAMPLES];
      if (flags & IPFLAG_FIND_FIRST_WM_PEAK) {
        int i;
        for (i = 0; i < numberOfSamples; i++) {
          if (i < numberOfSamples - 1 && i > 0)
            dm2[i] = dm[i + 1] - dm[i - 1];
          else
            dm2[i] = 0;
        }
        if (vno == Gdiag_no) {
          char fname[STRLEN];
          sprintf(fname, "v%d.%2.0f.dm.log", Gdiag_no, sigma * 100);
          fp = fopen(fname, "w");
          for (i = 0; i < numberOfSamples; i++) fprintf(fp, "%f %f\n", dm[i], dm2[i]);
          fclose(fp);
          fp = NULL;
          DiagBreak();
        }
      }

      if (max_mag_val > 0 && max_mri / 1.15 > max_mag_val) {

        float peak = 0.0f, outside = 1.0f;
      
        int i;
        for (i = 0; i < numberOfSamples; i++) {
          if (i == Gdiag_no2) DiagBreak();
          
          // Find a local maxima, where the slope changes from positive to zero or negative
          if (dm[i] > 0) continue;

          peak    = dm[i];
          outside = 0.0;
          int num = 0;
          
          int const lo = MAX(0, i - whalf);          
          int const hi = MIN(i + whalf + 1, numberOfSamples);
          int i1;
          for (i1 = lo; i1 < hi; i1++) {
            outside += dm[i1];
            num++;
            if (dm[i1] < peak) break;  // not a local maxima in the negative direction
                // If i1 <  i then dm[i1] is positive, and peak is negative, so this test is false
                // If i1 == i then dm[i] == dm[i1], so this test fails
                // It i1 >  i then this is searching for the first following slope that is even steeper down
          }
          
          outside /= num;

          if ((peak < 0) && (i1 > i + whalf))  // found a local maximum that is not a flat region of 0
            break;
        }

        if (i < numberOfSamples - whalf && peak / outside > 1.5)  // it was a local max - set the target to here
        {
          if (vno == Gdiag_no)
            printf("v %d: resetting target to local max at %2.2f: I=%d, peak=%2.2f, outside=%2.2f, ratio=%2.2f\n",
                   vno,
                   dists[i],
                   (int)mri[i],
                   peak,
                   outside,
                   peak / outside);
          max_mag_val  = mri[i];            // beware - this is finding the highest vno that gets here
          max_mag      = fabs(dm[i]);       // but getting here depends on the previously found max_mag_val !
          max_mag_dist = dists[i];          // so this code can not be parallelized as is...
        }
        
        else if (flags & IPFLAG_FIND_FIRST_WM_PEAK)  // not a local max in 1st derivative - try second */
        {
          for (i = 0; i < numberOfSamples; i++) {
            if (i == Gdiag_no2) DiagBreak();

            if (dm2[i] >= 0) continue;

            peak = dm2[i];
            int num = 0;
            outside = 0.0;
            
            int const lo = MAX(0, i - whalf);          
            int const hi = MIN(i + whalf + 1, numberOfSamples);
            int i1;
            for (i1 = lo; i1 < hi; i1++) {
              outside += dm2[i1];
              num++;
              if (dm2[i1] < dm2[i]) break;  // not a local maxima in the negative direction
            }
            outside /= num;
            double val          = mri[i];
            double next_val     = mri[i + 1];
            double previous_val = mri[i - 1];
            // make sure it is in feasible range
            if ((previous_val > inside_hi) || (previous_val < border_low) || (next_val > outside_hi) ||
                (next_val < outside_low) || (val > border_hi) || (val < border_low))
              continue;

            if ((peak < 0) && (i1 > i + whalf) && mri[i1])  // found a local maximum that is not a flat region of 0
              break;
          }

          if (i < numberOfSamples - whalf && peak / outside > 1.5)  // it was a local max - set the target to here
          {
            if (vno == Gdiag_no)
              printf(
                  "!!!!!!!!! v %d: resetting target to local max at in second derivative %2.2f: I=%d, peak=%2.2f, "
                  "outside=%2.2f, ratio=%2.2f\n",
                  vno,
                  dists[i],
                  (int)mri[i],
                  peak,
                  outside,
                  peak / outside);

            max_mag = dm[i];
            
            int i1;
            for (i1 = i + 1; i1 < numberOfSamples; i1++)  // search forward for largest (negative) derivative
              if (max_mag > dm[i1])           // previous one was largest negative one
                break;
            
            if (i1 < numberOfSamples) i = i1 - 1;
            
            max_mag_val  = mri[i];
            max_mag      = fabs(dm[i]);
            max_mag_dist = dists[i];
          }
        }
        
        if (vno == Gdiag_no) DiagBreak();
      }
    }

    if (which == GRAY_CSF && local_max_found == 0 && max_mag_dist > 0) {
      
      /* check to make sure it's not ringing near the gray white boundary,
         by seeing if there is uniform stuff outside that could be gray matter.
      */
      int allgray = 1;

      float outlen;
      for (outlen = max_mag_dist; outlen < max_mag_dist + 2; outlen += STEP_SIZE) {
        double const x = v->x + v->nx * outlen;
        double const y = v->y + v->ny * outlen;
        double const z = v->z + v->nz * outlen;
        double xw,yw,zw;
        MRIS_useRAS2VoxelMap(sras2v_map, mri_brain,x, y, z, &xw, &yw, &zw);
        double val;
        MRIsampleVolume(mri_brain, xw, yw, zw, &val);
        if ((val < outside_hi /*border_low*/) || (val > border_hi)) {
          allgray = 0;
          break;
        }
      }
      
      if (allgray) {
        if (Gdiag_no == vno)
          printf(
              "v %d: exterior gray matter detected, "
              "ignoring large gradient at %2.3f (I=%2.1f)\n",
              vno,
              max_mag_dist,
              max_mag_val);

        max_mag_val = -10; /* don't worry about largest gradient */
        max_mag_dist = 0;
        num_changed++;
      }
    }

    if (max_mag_val > 0) /* found the border value */
    {
      if (local_max_found) {
        ngrad_max++;
      }
      else {
        ngrad++;
      }
      if (max_mag_dist > 0) {
        nout++;
        nfound++;
        mean_out += max_mag_dist;
      }
      else {
        nin++;
        nfound++;
        mean_in += -max_mag_dist;
      }

      if (max_mag_val < border_low) {
        max_mag_val = border_low;
      }

      mean_dist += max_mag_dist;

      v->val  = max_mag_val;
      v->mean = max_mag;
      
      mean_border += max_mag_val;
      total_vertices++;
      
      v->d = max_mag_dist;
      v->marked = 1;
    }
    else /* couldn't find the border value */
    {
      if (min_val < 1000) {
        nmin++;
        v->d = min_val_dist;
#if 0
        if (min_val > border_hi)  /* found a low value, but not low enough */
        {
          min_val = border_hi ;
        }
        else if (min_val < border_low)
        {
          min_val = border_low ;
        }
#else
        if (min_val < border_low) {
          min_val = border_low;
        }
#endif
        v->val = min_val;
        mean_border += min_val;
        total_vertices++;
        v->marked = 1;
      }
      else {
        /* don't overwrite old target intensity if it was there */
        /*        v->val = -1.0f ;*/
        v->d = 0;
        if (v->val < 0) {
          nalways_missing++;
#if 0
          v->val = (border_low+border_hi)/2 ;
#endif
          v->marked = 0;
        }
        else {
          v->marked = 1;
        }
        nmissing++;
      }
    }
    if (vno == Gdiag_no)
      fprintf(stdout,
              "v %d, target value = %2.1f, mag = %2.1f, dist = %2.2f, %s\n",
              Gdiag_no,
              v->val,
              v->mean,
              v->d,
              local_max_found ? "local max" : max_mag_val > 0 ? "grad" : "min");
#if 0
    if (vno == 44289 || vno == 91080 || vno == 92286 || vno == 46922)
      fprintf(stdout, "v %d, target value = %2.1f, mag = %2.1f, dist=%2.2f\n",
              Gdiag_no, v->val, v->mean, v->d) ;
#endif
  
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  mean_dist   /= (float)(total_vertices - nmissing);
  mean_border /= (float)total_vertices;

  if (nin > 0) {
    mean_in   /= (float)nin;
  }
  if (nout > 0) {
    mean_out /= (float)nout;
  }

#if 0
  MRISsoapBubbleVals(mris, 100) ;
#endif

  /*  MRISaverageVals(mris, 3) ;*/
  fprintf(stdout,
          "mean border=%2.1f, %d (%d) missing vertices, mean dist %2.1f "
          "[%2.1f (%%%2.1f)->%2.1f (%%%2.1f))]\n",
          mean_border,
          nmissing,
          nalways_missing,
          mean_dist,
          mean_in,
          100.0f * (float)nin  / (float)nfound,
          mean_out,
          100.0f * (float)nout / (float)nfound);

  fprintf(stdout,
          "%%%2.0f local maxima, %%%2.0f large gradients "
          "and %%%2.0f min vals, %d gradients ignored\n",
          100.0f * (float)ngrad_max / (float)mris->nvertices,
          100.0f * (float)ngrad     / (float)mris->nvertices,
          100.0f * (float)nmin      / (float)mris->nvertices,
          num_changed);

  if (log_fp) {

    fprintf(log_fp,
            "mean border=%2.1f, %d (%d) missing vertices, mean dist %2.1f "
            "[%2.1f (%%%2.1f)->%2.1f (%%%2.1f))]\n",
            mean_border,
            nmissing,
            nalways_missing,
            mean_dist,
            mean_in,
            100.0f * (float)nin  / (float)nfound,
            mean_out,
            100.0f * (float)nout / (float)nfound);

    fprintf(log_fp,
            "%%%2.0f local maxima, %%%2.0f large gradients "
            "and %%%2.0f min vals, %d gradients ignored\n",
            100.0f * (float)ngrad_max / (float)mris->nvertices,
            100.0f * (float)ngrad     / (float)mris->nvertices,
            100.0f * (float)nmin      / (float)mris->nvertices,
            num_changed);
  }

  MRIS_freeRAS2VoxelMap(&sras2v_map);

  return (NO_ERROR);
}


