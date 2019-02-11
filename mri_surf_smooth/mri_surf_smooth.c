/**
 * @file  mri_surf_smooth.c
 * @brief smoothing along the surfaces with given neighborhood size and radial smoothing
 *
 *
 */
/*
 * Original Author: Anna Blazejewska
 * CVS Revision Info:
 *    $Author:  $
 *    $Date:  $
 *    $Revision:  $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "tags.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"
#include "fsenv.h"
#include "registerio.h"

// TODO time measuring, remove after testing is finished
#include <sys/time.h>
#include <time.h>
#include <unistd.h>


static int parse_commandline(int argc, char **argv);
static void print_help(void);
static int get_neighborhood(MRI_SURFACE *mris, int v_id, int size, int *neighbors, int *nb_weights,  int *current_nb_count);

int main(int argc, char *argv[]) ;

char *Progname, name[STRLEN], hemi[2] = "lh", text[STRLEN] = "", depth_fname[STRLEN], reg_fname[STRLEN];
int nb_size = 0, rad_size = 0,  surf_num = 11, vol_map = 1, nb_count = 0, nb_count_w = 0, rad_pve = 0, rad = 0, rad_nb = 0;
float rad_surf = 1.0;
FSENV *fsenv;
char *subject = NULL;
static char vcid[] = "$Id: mris_surf_smooth.c,v 1.30 2014/01/21 18:48:21 Anna Exp $";


int main(int argc, char *argv[]) {
  char out_fname[STRLEN], path[STRLEN], in_fname[STRLEN], map_fname[STRLEN];
  int  f, i = 0, j, t;
  clock_t begin;

/*
  char cmdline[CMD_LINE_LEN] ;

  make_cmd_version_string
  (argc, argv,
   "$Id: mris_smooth.c,v 1.30 2014/01/21 18:48:21 fischl Exp $",
   "$Name:  $", cmdline);
*/
  /* rkt: check for and handle version tag */
	/*
  nargs = handle_version_option
          (argc, argv,
           "$Id: mris_smooth.c,v 1.30 2014/01/21 18:48:21 fischl Exp $",
           "$Name:  $");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
*/

  fsenv = FSENVgetenv();
  parse_commandline(argc, argv);

///////////////////////////////////////////////////////////////////////////////
  MRI_SURFACE *mris[surf_num];
  MRI *map[surf_num], *pve[surf_num], *output, *debug;

  sprintf(path, "%s/%s", fsenv->SUBJECTS_DIR, subject);
  printf("Processing subject: %s\n", path);

  for (f=0; f<surf_num; f++) {
    // read surface data
    sprintf(in_fname, "%s/recon/surf/%s.midgray.%2.2d", path, hemi, f) ;

    printf("Reading: %s\n", in_fname);
    mris[f] = MRISfastRead(in_fname);
    if (!mris[f]) ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",  Progname, in_fname);
  }

  if (mris[5]->nvertices >500000)
	  ErrorExit(ERROR_BADPARM, "Number vertices %d is higher than 500000",  Progname, in_fname);

  // TODO parameter!
  //const long int u = mris[5]->nvertices;
  static int neighbors[500000][350];
  static int nb_weights[500000][350];
  static int current_count[500000];

  begin = clock();

  // neighborhood option
  if (nb_size>0) {
    for (i=0; i<mris[5]->nvertices; i++) {
      for (j=0; j<mris[5]->nvertices; j++) mris[5]->vertices[j].marked = 0;
      //this is HACK: size+1 so the 1st neighbors are marked with the nb size
      get_neighborhood(mris[5], i, nb_size+1, neighbors[i], nb_weights[i], &current_count[i]);
      for (t=0; t<current_count[i]; t++) {
    	nb_weights[i][t] = mris[5]->vertices[neighbors[i][t]].marked;
    	// re-mark the neighbour to the right weight
     	if (neighbors[i][t] == i) nb_weights[i][t] = nb_size;
      }
    }
  }

  printf("time %f\n", (clock()-begin)/(float)CLOCKS_PER_SEC);

//-------------------------------------------------------------------------------------
//// voxels within NB WEIGHTED or NOT counting only & exit ----------------------------
//-------------------------------------------------------------------------------------
	if (nb_count == 1 || nb_count_w == 1) {
		long int nb_index[350];
		float unique = 0;
		float nb_w[350];

		// read indexed overlay data
		for (f=0; f<surf_num; f++) {
			sprintf(map_fname, "%s/map_surf_index/%s.%s.midgray_index.%2.2d.mgz", path, hemi, name, f) ;
			printf("Reading: %s\n", map_fname);
			map[f] = MRIread(map_fname);
			if (map[f]==NULL) ErrorExit(ERROR_NOFILE, "%s: could not read mri file %s",  Progname, map_fname);
		}
		sprintf(out_fname, "%s/map_surf_index/%s_%s_%s_midgray_index_%s_voxn_wt.csv",  path, subject, hemi, name, text) ;
		printf("Writting to: %s\n", out_fname);

		FILE *pFile;
		pFile = fopen (out_fname,"w");
		for (f=0; f<surf_num; f++) {
			unique = 0;
			for (i=0; i<mris[f]->nvertices; i++) {
				// for each neighbor found voxel id
				for (j=0; j<current_count[i]; j++) {
					nb_index[j] = (long int)MRIgetVoxVal(map[f], neighbors[i][j], 0, 0, 0);
					nb_w[j] = 1.0;
					if (nb_count_w == 1) nb_w[j] = nb_w[j]/(nb_size - nb_weights[i][j]+1);
					//printf("%f ", nb_w[j]);
				}

				// count the unique ones
				for (j=0; j<current_count[i]; j++) {
					unique =  unique + nb_w[j];
					for (t=0; t<j; t++) {
						if (nb_index[j] == nb_index[t]) {
							unique = unique - nb_w[j];
							break;
						}
					}
				}
				// cleaning
				for (j=0; j<current_count[i]; j++) {
					nb_index[j] = 0;
					nb_w[j] = 0.0;
				}
			}
			printf("%f ", (float)unique/mris[f]->nvertices);
			fprintf (pFile, "%f;", (float)unique/mris[f]->nvertices);
		}
		printf("\n");
		fclose (pFile);
		//cleaning
		for (f=0; f<surf_num; f++) MRISfree(&mris[f]);
		// TIME
		printf("time %f\n", (clock()-begin)/(float)CLOCKS_PER_SEC);
		return(NO_ERROR);
	}
//-------------------------------------------------------------------------------------
//// ----------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------

  // read overlay data
  for (f=0; f<surf_num; f++) {
  	sprintf(map_fname, "%s/map_surf/%s.%s.midgray.%2.2d.*.mgz", path, hemi, name, f) ;
    printf("Reading: %s\n", map_fname);
    map[f] = MRIread(map_fname);
    if (map[f]==NULL) ErrorExit(ERROR_NOFILE, "%s: could not read mri file %s",  Progname, map_fname);
  }
  printf("%d %d %d %d\n", map[0]->width, map[0]->height, map[0]->depth, map[0]->nframes);


  // if radial partial volume smoothing, read partial volume overlays
  if (rad_pve == 1 || rad == 1) {
    for (f=0; f<surf_num; f++) {
	  sprintf(map_fname, "%s/partial_volume/%s.%s.midgray_pv.%2.2d.mgz", path, hemi, name, f) ;
	  printf("Reading: %s\n", map_fname);
	  pve[f] = MRIread(map_fname);
	  if (pve[f]==NULL) ErrorExit(ERROR_NOFILE, "%s: could not read mri file %s",  Progname, map_fname);
    }
  }


  // for each surface
  for (f=0; f<surf_num; f++) {
    float count = 0, count_zero = 0;
    float val[map[0]->nframes];

    // initializing the output
	  output = MRIclone(map[f], NULL);
	  debug = MRIclone(map[f], NULL);

	  // if radial smoothing requested and not possible
	  if ((rad_size>0) && ((f-rad_size)<0 || (f+rad_size)>=surf_num)) printf("Not enough layers to perform radial smoothing\n");

	  // for each vertex of the surface
	  for (i=0; i<mris[f]->nvertices; i++) {

	  	// initialize
			for (t=0; t<map[f]->nframes; t++) val[t] = 0;
			count = 0;
			count_zero = 0;

			// RAD / RAD-PVE --- same result for each surface -> 1 output file -------------------
 	    if (rad == 1 || rad_pve == 1) {
 	    	float w = 1.0;
			  for (j=0; j<surf_num; j++) {
			  	// detecting vertices = 0 on all surfaces
			  	// for CSF
			  	//if (MRIgetVoxVal(pve[j], i, 0, 0, 0) == 1) count_zero++;
			  	if (MRIgetVoxVal(pve[j], i, 0, 0, 0) == 0) count_zero++;
			  	w = 1.0;
			  	// W for CSF
				  //if (rad_pve == 1) w = 1-MRIgetVoxVal(pve[j], i, 0, 0, 0);
				  // W for GM
				  if (rad_pve == 1) w = MRIgetVoxVal(pve[j], i, 0, 0, 0);
			  	for (t=0; t<map[f]->nframes; t++) val[t] += MRIgetVoxVal(map[j], i, 0, 0, t)*w;
			  	count += w;
  		  }
			  if (count_zero == surf_num) {
			  	for (t=0; t<map[f]->nframes; t++) val[t] = 0;
			  	count = 1;
			  }
			// NB-SIZE -------------------------------------------------------------------------
 	    } else if (nb_size>0) {
				// calculate smoothed value using marked neighbors
				for (j=0; j<current_count[i]; j++) {
					for (t=0; t<map[f]->nframes; t++)
						val[t] += MRIgetVoxVal(map[f], neighbors[i][j], 0, 0, t)/(nb_size-nb_weights[i][j]+1);
						//val += ((float *)map[f]->slices[t][0])[j]/(nb_size-mris[f]->vertices[neighbors[j]].marked+1); //pointer version
					count += 1.0/(nb_size - nb_weights[i][j]+1);
				}
			}
			// RAD-SIZE -------------------------------------------------------------------------
			if (rad_nb == 1) {
			  for (j=0; j<surf_num; j++) {
					for (t=0; t<map[f]->nframes; t++)
						val[t] += MRIgetVoxVal(map[f-j], i, 0, 0, t);
				count++;
				}
			}
 	    /*
			// RAD-SIZE -------------------------------------------------------------------------
			if ((rad_size>0) && ((f-rad_size)>=0 && (f+rad_size)<surf_num)) {
				for (j=1; j<=rad_size; j++) {
					for (t=0; t<map[f]->nframes; t++)
					val[t] += MRIgetVoxVal(map[f+j], i, 0, 0, t)/j + MRIgetVoxVal(map[f-j], i, 0, 0, t)/j;
				count += 2/j;
				}
			}


			// RAD-SURF only if surf_num>0 ---------------------------------------------------------
			for (j=0; j<surf_num; j++) {
				for (t=0; t<map[f]->nframes; t++)
				//skip negative ones
				if ((1-fabs(j-f)*rad_surf)>0) { printf("im doing stuff"); val[t] += MRIgetVoxVal(map[j], i, 0, 0, t)*(1-fabs(j-f)*rad_surf);}
				if ((1-fabs(j-f)*rad_surf)>0) count += (1-fabs(j-f)*rad_surf);
			}
			 */

			// assign the new value to the output overlay
			for (t=0; t<map[f]->nframes; t++) {
				//printf("%f %f | ", val[t], val[t]/count);
				MRIsetVoxVal(output,i,0,0,t, val[t]/count);
			}
			MRIsetVoxVal(debug,i,0,0,0, count);
		//MRIdbl2ptr(val[t]/count, output->slices[t][0]+i, output->type); //pointer version
	  }



    // TODO: move it to the right place!
    if (rad_nb == 1) sprintf(text, "dyad_nb%d", nb_size);

	  // write output
	  sprintf(out_fname, "%s/map_surf_smooth/%s.%s.midgray_smooth_%s.%2.2d.mgz", path, hemi, name, text, f);
	  printf("Writing: %s\n", out_fname);
	  MRIwrite(output, out_fname);

	  //write debug
/*
	  sprintf(out_fname, "%s/map_surf_c/%s.%s.midgray_smooth_%s.%2.2d_W_csf_test.mgz", path, hemi, name, text, f);
	  printf("Writing: %s\n", out_fname);
	  MRIwrite(debug, out_fname);
*/
	  // freeing memory
	  MRIfree(&output);

	  // TIME
	  printf("time %f\n", (clock()-begin)/(float)CLOCKS_PER_SEC);
	  if (rad == 1 || rad_pve == 1) return(NO_ERROR);
  }

////////////////////////////////////////////////////////////////////////////////
// MAPPING TO THE VOLUME
// only if the depth map and the registration matrix has been specified
////////////////////////////////////////////////////////////////////////////////
  if (vol_map) {
    MRI *depth, *orig;

    //char *volregfile = NULL;
    float ipr, bpr, intensity;
    int float2int, err;
    MATRIX *Ma2vTKR;

    // read original data
    sprintf(in_fname, "%s/recon/mri/orig.mgz", path);
		printf("Reading: %s\n", in_fname);
		orig = MRIread(in_fname);
		if (orig==NULL) ErrorExit(ERROR_NOFILE, "%s: could not read mri file %s", Progname, in_fname);

		// read depth data
    sprintf(map_fname, "%s/data/%s", path, depth_fname );
		printf("Reading: %s\n", map_fname);
		depth = MRIread(map_fname);
		if (depth==NULL) ErrorExit(ERROR_NOFILE, "%s: could not read mri file %s", Progname, map_fname);

		// read in the tkregister registration
    sprintf(in_fname, "%s/data/%s", path, reg_fname);
		err = regio_read_register(in_fname, &subject, &ipr, &bpr, &intensity, &Ma2vTKR, &float2int);
		if (err) ErrorExit(ERROR_NOFILE, "%s: could not read file %s", Progname, in_fname);

		/////////////////////////////////////////////////////////////////////////////
  }

  if (rad_pve == 1)	for (f=0; f<surf_num; f++) MRIfree(&pve[f]);
  // freeing even more memory
  for (f=0; f<surf_num; f++) {
		MRIfree(&map[f]);
		MRISfree(&mris[f]);
  }

  return(NO_ERROR);
}


// gives the list of neighbors of the certain vertex, for a given neighborhood size
static int get_neighborhood(MRI_SURFACE *mris, int v_id, int size, int *neighbors, int *nb_weights, int *current_nb_count) {
  int j;
  // spread
  if (size>0) {
		//if (size==4)
		for (j=0; j<mris->vertices[v_id].vnum; j++ ) {
			//printf("NUM = %d ", mris->vertices[v_id].vnum);
			get_neighborhood(mris, mris->vertices[v_id].v[j], size-1, neighbors, nb_weights, current_nb_count);
		}
		//printf("M=%d ", mris->vertices[v_id].marked);
		// for new ones
		if (mris->vertices[v_id].marked == 0) {
			neighbors[*current_nb_count] = v_id;
			(*current_nb_count)++;
			//mris->vertices[v_id].marked = 1;
		}
		if (mris->vertices[v_id].marked < size) {
			mris->vertices[v_id].marked = size;
		}
  }
  return(NO_ERROR);
}

static int parse_commandline(int argc, char **argv) {
  int nargc , nargsused;
  char **pargv, *option ;
  printf("%d\n", argc);
  if (argc <= 1) print_help();
  nargc = argc;
  pargv = argv;
  while (nargc > 0) {
    option = pargv[0];
    nargc -= 1;
    pargv += 1;
    nargsused = 0;
    if (!stricmp(option, "-nb")) {
      if (nargc < 1) ErrorExit(ERROR_BADPARM, "Flag %s needs an argument\n", option);
      nb_size = atoi(pargv[0]) ;
      nargsused = 1;
	    sprintf(text, "nb%d", nb_size);
      if (nb_size > 10) ErrorExit(ERROR_BADPARM, "Neighborhood size = %d exceeds the limit of 10\n", nb_size);
      else fprintf(stderr, "Using neighborhood size = %d\n", nb_size);
    } else if (!stricmp(option, "-nb-count")) {
      nb_count = 1;
  	  fprintf(stderr, "Counting avg number of neighbors for each layer\n");
  	  sprintf(text, "voxn");
    } else if (!stricmp(option, "-nb-count-w")) {
      nb_count_w = 1;
  	  fprintf(stderr, "Counting avg number of neighbors for each layer weighted by their contribution\n");
  	  sprintf(text, "voxn_wt");
    } else if (!stricmp(option, "-rad-nowt")) {
      rad = 1;
      // label for the output file name
      sprintf(text, "rad");
      nargsused = 1;
   	  fprintf(stderr, "Radial smoothing with no weights\n");
    } else if (!stricmp(option, "-rad-nb")) {
      rad_nb = 1;
      // label for the output file name
      sprintf(text, "rad");
      nargsused = 1;
   	  fprintf(stderr, "Radial smoothing with no weights\n");
    } else if (!stricmp(option, "-rad-size")) {
      if (nargc < 1) ErrorExit(ERROR_BADPARM, "Flag %s needs an argument\n", option);
      rad_size = atoi(pargv[0]);
      // label for the output file name
	    sprintf(text, "rad%d", rad_size);
      nargsused = 1;
    } else if (!stricmp(option, "-rad-surf")) {
      if (nargc < 1) ErrorExit(ERROR_BADPARM, "Flag %s needs an argument\n", option);
      if (rad_surf < 0 || rad_surf > 1) ErrorExit(ERROR_BADPARM, "Radial smoothing slope = %f exceeds the range of [0,1]\n", rad_surf);
      else {
        rad_surf = atof(pargv[0]);
    	fprintf(stderr, "Using radial smoothing slope = %f\n", rad_surf);
      }
      // label for the output file name
  	  sprintf(text, "rad_slope%1.2f", rad_surf);
      nargsused = 1;
    } else if (!stricmp(option, "-rad-pve")) {
      rad_pve = 1;
   	  fprintf(stderr, "Radial smoothing with partial volume weights\n");
  	  sprintf(text, "rad-pve");
//--------- obligatory ------------------------------------------------------------
    } else if (!stricmp(option, "-in")) {
      if (nargc < 1) ErrorExit(ERROR_BADPARM, "Flag %s needs an argument\n", option);
  	  strcpy(name, pargv[0]);
      nargsused = 1;
  	  fprintf(stderr, "Name: %s\n", name);
    } else if (!stricmp(option, "-h")) {
      if (nargc < 1) ErrorExit(ERROR_BADPARM, "Flag %s needs an argument\n", option);
  	  strcpy(hemi, pargv[0]);
      nargsused = 1;
  	  if (strcmp(hemi,"lh") && strcmp(hemi,"rh")) ErrorExit(ERROR_BADPARM, "Hemi = %s, must be lh or rh\n",hemi);
  	  fprintf(stderr, "Hemisphere: %s\n", hemi);
    } else if (!stricmp(option, "-s")){
      if (nargc < 1) ErrorExit(ERROR_BADPARM, "Flag %s needs an argument\n", option);
      subject = pargv[0];
      nargsused = 1;
    } else if (!stricmp(option, "-n")) {
      if (nargc < 1) ErrorExit(ERROR_BADPARM, "Flag %s needs an argument\n", option);
      surf_num = atoi(pargv[0]) ;
      nargsused = 1;
      if (surf_num > 11) ErrorExit(ERROR_BADPARM, "Number of layers %d exceeds the limit of 11\n", surf_num);
      else fprintf(stderr, "Number of layers = %d\n", surf_num);
    } else if (!stricmp(option, "-depth")) {
      if (nargc < 1) ErrorExit(ERROR_BADPARM, "Flag %s needs an argument\n", option);
	  strcpy(depth_fname, pargv[0]);
      nargsused = 1;
	  fprintf(stderr, "Depth file name: %s\n", depth_fname);
    } else if (!stricmp(option, "-reg")) {
      if (nargc < 1) ErrorExit(ERROR_BADPARM, "Flag %s needs an argument\n", option);
	  strcpy(reg_fname, pargv[0]);
	  nargsused = 1;
	  fprintf(stderr, "Registration file name: %s\n", reg_fname);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }

  // validate dependent parameters
  if ((nb_count == 1) && (nb_count_w == 1)) ErrorExit(ERROR_BADPARM, "NB count can be weighted or not waited, but not both in the same time..\n");
  if ((rad == 1) && ((rad_surf < 1)||(rad_size > 0)||(rad_pve == 1))) ErrorExit(ERROR_BADPARM, "Flag -rad cannot be used in combination with any other -rad-x option \n");
  if ((rad_size > 0) && (rad_surf < 1)) ErrorExit(ERROR_BADPARM, "Flag -rad-surf cannot be used in combination with flag -rad-size\n");
  if ((rad_pve == 1) && ((rad_size > 0)||(rad_surf < 1))) ErrorExit(ERROR_BADPARM, "Flag -rad-pve cannot be used in combination with any other -rad-x option\n");
  if (rad_size > ((surf_num-1)/2)) ErrorExit(ERROR_BADPARM, "Radial smoothing size = %d exceeds the limit determined by the number of surfaces %d\n", rad_size, surf_num);
  else fprintf(stderr, "Using radial smoothing size = %d\n", rad_size);
  if (!(strlen(depth_fname) != 0 && strlen(reg_fname) != 0)) {
	vol_map = 0;
    fprintf(stderr, "Volume mapping will not be performed\n");
  }
  if ((nb_count == 1) && (nb_size <= 0))
    fprintf(stderr, "Neighborhood size has to be specified to count voxels laying within this neighborhood\n");
  return(0);
}


static void print_help(void) {
  printf("\n%s\n\n",vcid);
  printf(
  "Smoothes the EPI data mapped to the surfaces with given tangmental and radial distance.\n"

  "mri_surf_smooth in_file map_file out_file -in NAME -s subject_dir [-nb-count -nb size -n number_of_surfaces -h hemisphere -depth depth_map -reg register.dat]\n"
    "-in NAME: base name of the input EPI file mapped to surfaces, see below for details\n"
    "-s subject_dir: subject directory\n"
    "-nb-count: ONLY count's (unique) voxels of EPI data which are included in given neighborhood for each surface\n"
  		"   requires specified nb size and folder sub/data/map_surf_index with voxel indexed surfaces\n"
    "-nb-count-w: ONLY count's (unique) voxels of EPI data which are included in given neighborhood for each surface\n"
  		"   weighted by their contribution\n"
  		"   requires specified nb size and folder sub/data/map_surf_index with voxel indexed surfaces\n"
  	"-nb size: neighborhood size [default: 0 = no smoothing]\n"
  		"   maximum ONE of two options below can be used: -rad-size or -rad-surf\n"
  	"-rad-size size: number of surfaces for radial smoothing (size*2+1)\n"
  		"   if a particular surface does not have enough neighboring ones,"
  		"   smoothing will skipped for this surface [default: 0 = smoothing off]\n"
    "-rad-surf slope: [0,1] slope of weighted smoothing across all the surfaces "
  		"[default: 1 = no smoothing]\n"
    "-rad-pve: radial smoothing across all the surfaces weighted by partial volume effect\n"
  		"   requires directory data/partial_volme with PV surface projections\n"
    "-rad-nowt: simple radial smoothing across all the surfaces with no weights\n"
  	"-n numbr_of_surfaces: [default: 11]\n"
    "-h lh|rh: [default: lh]\n\n"
    "-rad-nb: radial and nb\n"

  "Two files specified below should be in subjects 'data' directory:\n"
    "-depth depth_map: depth map for the hemisphere; if depth map is specified mapping to the volume will be performed using registration file specified with -reg flag\n"
    "-reg register.dat: registration file EPI -> structural used for mapping smooth surfaces back to the volume\n\n"

  "Assumed directory structure: subject/recon/surf, subject/map_surf \n"
    "Assumed file names: \n"
    "   mid surfaces (in surf directory): (l|r)h.midgray.xx.mgh \n"
    "   EPI mapped to the surfaces (in map_surf directory): (l|r)h.NAME.midgray.xx.mgh \n"
  	"\n");

  exit(1) ;
}

