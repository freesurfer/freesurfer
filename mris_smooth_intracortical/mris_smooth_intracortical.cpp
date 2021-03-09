/**
 * @brief smoothing along the cortical surface meshes with a given tangential (neighborhood size) and radial (number of surface meshes) extent
 *
 * This method has been described in:
 * Blazejewska, AI, Fischl, B, Waldab, LL, Polimeni
 * Intracortical smoothing of small-voxel fMRI data can provide increased detection power without spatial resolution losses compared to conventional large-voxel fMRI data
 * 2019, Neuroimage 189 (601-614)
 * https://doi.org/10.1016/j.neuroimage.2019.01.054
 */

/*
 * Original Author: Anna I. Blazejewska
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

/*!
 \file mri_surf_smooth.cpp
 \brief Smoothing along the cortical surface meshes with a given tangential and radial extent ((neighborhood size &number of surface meshes).
 \author Anna I. Blazejewska
 */

/*
 BEGINHELP
 TODO
 ENDHELP
 */

/*
 BEGINUSAGE
 TODO
 ENDUSAGE
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <glob.h>
#include <libgen.h>

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

#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#define MAX_NB (6) // the reasonable number of neighbors (tan size)
#define MAX_SURF (20) // max number of surfaces
#define MAX_VERTICES (1000000) // max number of vertices
#define SEP "/"

#ifndef GLOB_PERIOD
#define GLOB_PERIOD 0
#endif

int main(int argc, char *argv[]);
static void calculate_nb_weights(float *nb_weights, int nb_num, int *hops);

static int parse_commandline(int argc, char **argv);
static void print_help(void);
static void check_options(void);
static void print_usage(void);
static void print_version(void);



const char *Progname = NULL;

char surf_path[STRLEN], over_path[STRLEN], out_path[STRLEN], surf_name[STRLEN], over_name[STRLEN], surf_dir[STRLEN], over_dir[STRLEN], out_dir[STRLEN], out_name[STRLEN];
int surf_num = 0, over_num = 0, nb_rad = 0, ic_size = 1, ic_start = 0, nb_wf = 0; //nb_wf = 0 (gauss)

int main(int argc, char *argv[]) {
	Progname = argv[0];

	int f, v, n, t;
	clock_t begin;
	begin = clock();

	//fsenv = FSENVgetenv();
	parse_commandline(argc, argv);

	// get input files paths
	glob_t over_list, surf_list;
	glob(strcat(strcat(strcpy(over_path, over_dir), SEP), over_name), GLOB_PERIOD, NULL, &over_list);
	glob(strcat(strcat(strcpy(surf_path, surf_dir), SEP), surf_name), GLOB_PERIOD, NULL, &surf_list);
	surf_num = surf_list.gl_pathc; over_num = over_list.gl_pathc;

	// check consistency of other options with number of files and exits with error if any inconsistency is detected
	check_options();

	MRI_SURFACE *surf[ic_size];
	MRI *over[ic_size], *output[ic_size];

	// read only surfaces/overlays that are included in the radial extent of the kernel
	for (f = ic_start; f < ic_size+ic_start; f++) {
		// read overlays
		printf("Reading: %s\n", over_list.gl_pathv[f]);
		over[f-ic_start] = MRIread(over_list.gl_pathv[f]);
		if (over[f-ic_start] == NULL) ErrorExit(ERROR_NOFILE, "%s: could not read MRI overlay file %s.\n", Progname,	over_list.gl_pathv[f]);
		if ((f>ic_start) && (over[f-ic_start]->nframes != over[0]->nframes))
			ErrorExit(ERROR_BADPARM, "Number of frames in overlay %s = %d is different than number of frames in overlay %s = %d.\n", over_list.gl_pathv[f], over[f-ic_start]->nframes, over_list.gl_pathv[0], over[0]->nframes);
		// read surfaces, currently only required for tangential smoothing.. this will change in the future
		if (nb_rad) {
			printf("Reading: %s\n", surf_list.gl_pathv[f]);
			surf[f-ic_start] = MRISread(surf_list.gl_pathv[f]);
			if (!surf[f-ic_start]) ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s.\n", Progname, surf_list.gl_pathv[f]);
			// all surfaces must have the same number of vertices smaller than MAX_VERTICES
			if (surf[f-ic_start]->nvertices > MAX_VERTICES) ErrorExit(ERROR_BADPARM, "Number of vertices %d is higher than %d.\n", surf_list.gl_pathv[0], MAX_VERTICES);
			if ((f>ic_start) && (surf[f-ic_start]->nvertices != surf[0]->nvertices))
				ErrorExit(ERROR_BADPARM, "Number of vertices in surface %s = %d is different than number of vertices in surface %s = %d.\n", surf_list.gl_pathv[f], surf[f-ic_start]->nvertices, surf_list.gl_pathv[0], surf[0]->nvertices);
			if (over[f-ic_start]->width != surf[f-ic_start]->nvertices)
				ErrorExit(ERROR_BADPARM, "Number of data points in overlay %s = %d is different than number of vertices in surface %s = %d.\n", over_list.gl_pathv[f], over[f-ic_start]->width, surf_list.gl_pathv[f], surf[f-ic_start]->nvertices);
		}
	}


	// if no output dir/name given - set based on 1st overlay
	if (strlen(out_dir) == 0) strcpy(out_dir, over_dir);
	if (strlen(out_name) == 0) {
		char *dot, *name;
		name = basename(over_list.gl_pathv[0]);
		dot = strrchr(name, '.') ;
		strncpy(out_name, name, dot - name);
		sprintf(out_name, "%s.nb%d_rad%d-%d.mgz", out_name, nb_rad, ic_start, (ic_start+ic_size-1));
	}

	globfree(&surf_list);	globfree(&over_list);

	// tangential smoothing part ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if (nb_rad > 0) {
		static int neighbors[MAX_NEIGHBORS];
		static int hops[MAX_NEIGHBORS];
		static float nb_weights[MAX_NEIGHBORS];
		static int nb_num;

		// initializing the output with 0s
		for (f = 0; f < ic_size; f++) output[f] = MRIclone(over[f], NULL);
		float val, count;
		// for each surface & each vertex
		for (f = 0; f < ic_size; f++) {
			for (v = 0; v < surf[f]->nvertices; v++) {
				// neighborhood
				nb_num = MRISfindNeighborsAtVertex(surf[f], v, nb_rad, MAX_NEIGHBORS, neighbors, hops);
				//TODO more weights
				calculate_nb_weights(nb_weights, nb_num, hops);
				// for each frame in the overlay corresponding to this surface
				for (t = 0; t<over[f]->nframes; t++) {
					// for center vertex always weight 1.0
					val = MRIgetVoxVal(over[f], v,  0, 0, t);
					count = 1.0;
					for (n = 0; n < nb_num; n++) {
						val += (MRIgetVoxVal(over[f], neighbors[n],  0, 0, t)*nb_weights[n]);
						count += nb_weights[n];
					}
					MRIsetVoxVal(output[f],v,0,0,t, val/count);
				}
			}
		}
	}

	// radial smoothing part ///////////////////////////////////////////////////////////////////////////////////////
	// TODO in 2.0 add weights
	if (ic_size > 1) {
		// initializing outputs with input overlays
		if (nb_rad == 0) for (f = 0; f < ic_size; f++) output[f] = MRIcopy(over[f], NULL);

		float val;
		// only for ic smoothing so all overlays have the same number of frames
		for (t=0; t<over[0]->nframes; t++) {
			for (v=0; v<over[0]->width; v++) {
				val = 0.0;
				for (f = 0; f < ic_size;f++)
					val += MRIgetVoxVal(output[f], v,  0, 0, t);
				// write to output 0
				MRIsetVoxVal(output[0],v,0,0,t, val/ic_size);
			}
		}
	}

	// write an output overlay ///////////////////////////////////////////////////////////////////////////////////////
	sprintf(out_path, "%s%s%s", out_dir, SEP, out_name);
	printf("Saving result: %s\n", out_path);
	MRIwrite(output[0], out_path);

	// free
	for (f = 0; f < ic_size; f++) {
		MRISfree(&surf[f]);
		MRIfree(&over[f]);
		MRIfree(&output[f]);
	}

	printf("Finished in %.2f seconds\n", (double)(clock() - begin) / CLOCKS_PER_SEC);
	return (NO_ERROR);
}




static void calculate_nb_weights(float *nb_weights, int nb_num, int *hops) {
	int n;
	// gauss
	if (nb_wf == 0) {
		float sigma = nb_rad/2.3548;
		for (n = 0; n < nb_num; n++) nb_weights[n] = (1/(sigma*sqrt(2*PI)))*exp(-(hops[n]*hops[n])/(2*sigma*sigma));
		//fmax(nb_weights);
	// 1/NB
	} else if (nb_wf == 1)	{
			for (n = 0; n < nb_num; n++) {
			if (hops[n] == 0) nb_weights[n] = 1.0; // for the center vertex
			else nb_weights[n] = (float)(1.0/hops[n]);
		}
	}
	return;
}




/*!
\fn int parse_commandline(int argc, char **argv)
\brief Parses the command-line arguments
\param argc - number of command line arguments
\param argv - pointer to a character pointer
*/
static int parse_commandline(int argc, char **argv) {
	int nargc, nargsused;
	char **pargv, *option;
	if (argc <= 9) print_help();
	nargc = argc-1;
	pargv = argv+1;
	while (nargc > 0) {
		option = pargv[0];
		nargc -= 1;
		pargv += 1;
		nargsused = 0;

		if (!stricmp(option, "--surf_dir")) {
			if (nargc < 1) ErrorExit(ERROR_BADPARM, "Flag %s needs an argument\n", option);
			strcpy(surf_dir, pargv[0]);
			nargsused = 1;
		} else if (!stricmp(option, "--surf_name")) {
			if (nargc < 1) ErrorExit(ERROR_BADPARM, "Flag %s needs an argument\n", option);
			strcpy(surf_name, pargv[0]);
			nargsused = 1;
		} else if (!stricmp(option, "--overlay_dir")) {
			if (nargc < 1) ErrorExit(ERROR_BADPARM, "Flag %s needs an argument\n", option);
			strcpy(over_dir, pargv[0]);
			nargsused = 1;
		} else if (!stricmp(option, "--overlay_name")) {
			if (nargc < 1) ErrorExit(ERROR_BADPARM, "Flag %s needs an argument\n", option);
			strcpy(over_name, pargv[0]);
			nargsused = 1;
		} else if (!stricmp(option, "--output_dir")) {
			if (nargc < 1) ErrorExit(ERROR_BADPARM, "Flag %s needs an argument\n", option);
			strcpy(out_dir, pargv[0]);
			nargsused = 1;
		} else if (!stricmp(option, "--output_name")) {
			if (nargc < 1) ErrorExit(ERROR_BADPARM, "Flag %s needs an argument\n", option);
			strcpy(out_name, pargv[0]);
			nargsused = 1;
		} else if (!stricmp(option, "--tan-size")) {
			if (nargc < 1) ErrorExit(ERROR_BADPARM, "Flag %s needs an argument\n", option);
			nb_rad = atoi(pargv[0]);
			nargsused = 1;
		} else if (!stricmp(option, "--rad-size")) {
			if (nargc < 1) ErrorExit(ERROR_BADPARM, "Flag %s needs an argument\n", option);
			ic_size = atoi(pargv[0]);
			nargsused = 1;
		} else if (!stricmp(option, "--rad-start")) {
			if (nargc < 1) ErrorExit(ERROR_BADPARM, "Flag %s needs an argument\n", option);
			ic_start = atoi(pargv[0]);
			nargsused = 1;
		// nb weights
		} else if (!stricmp(option, "--tan-weights")) {
			if (nargc < 1) ErrorExit(ERROR_BADPARM, "Flag %s needs an argument\n", option);
			if (!stricmp(pargv[0], "gauss")) nb_wf = 0;
			else if (!stricmp(pargv[0], "distance")) nb_wf = 1;
			else fprintf(stderr, "Unknown value %s for flag %s, a default gaussian weighting function will be applied instead.\n", pargv[0], option);
			nargsused = 1;
		} else if (!stricmp(option, "--help")) {
			print_help();
			exit(0);
		} else fprintf(stderr, "Unrecognized flag/argument %s will be ignored.\n", option);
		nargc -= nargsused;
		pargv += nargsused;
	}
	return (0);
}

/*!
\fn static void check_options(void)
\brief Checks command-line options
*/
static void check_options(void) {
	// check if ranges of the parameters are correct
	if (surf_num < 1) ErrorExit(ERROR_BADPARM, "At least 1 input surface is required. \n");
	if (over_num < 1) ErrorExit(ERROR_BADPARM, "At least 1 input overlay is required. \n");
	if ((nb_rad < 0) || (nb_rad>MAX_NB)) ErrorExit(ERROR_BADPARM, "Tangential neighborhood radius = %d exceeds the range of [0, %d].\n", nb_rad, MAX_NB);
	if ((ic_size < 1) || (ic_size>surf_num)) ErrorExit(ERROR_BADPARM, "Radial smoothing extent = %d exceeds the range of [1,number of input surfaces = %d].\n", ic_size, surf_num);
	if ((ic_start < 0) || (ic_start>(surf_num-ic_size))) ErrorExit(ERROR_BADPARM, "Radial smoothing starting surface = %d exceeds the range of [0, (number of input surfaces - radial smoothing extent) = %d].\n", ic_start, surf_num-ic_size);

	// no smoothing case, default parameters
	if ((ic_size == 1) && (nb_rad == 0)) {
		fprintf(stderr, "Neighborhood radius = %d & radial extent = %d => no smoothing will be performed.\n", nb_rad, ic_size);
		exit(0);
	// tangential smoothing: one overlay & one surface
	} else if ((ic_size == 1) && (nb_rad > 0)) {
		if (surf_num != 1) ErrorExit(ERROR_BADPARM, "Number of input surfaces = %d. Tangential smoothing requires 1 input surface.\n", surf_num);
		if (over_num != 1) ErrorExit(ERROR_BADPARM, "Number of input overlays = %d. Tangential smoothing requires 1 input overlay.\n", over_num);
		fprintf(stderr, "Neighborhood radius = %d & radial extent = %d => tangential smoothing will be performed for an input overlay.\n", nb_rad, ic_size);
	// intracortical smoothing
	} else if ((ic_size > 1) && (nb_rad >= 0)) {
		if (surf_num != over_num) ErrorExit(ERROR_BADPARM, "Number of input surfaces = %d has to be the same as number of input overlays = %d for intracortical smoothing.\n", surf_num, over_num);
		fprintf(stderr, "Neighborhood radius = %d & radial extent = %d from a starting surface = %d => intracortical smoothing will be performed.\n", nb_rad, ic_size, ic_start);
	}


	return;
}

/*!
\fn static void print_usage(void)
\brief Prints usage and exits
*/
static void print_usage(void) {
  printf("USAGE: %s --surf_dir surfdir --surf_name surfname --overlay_dir overdir --overlay_name overname [--output_dir outdir --output_name outname --tan-size tansize --rad-size radsize --rad-start radstart]\n",Progname);
  printf("\n"
					"  --surf_dir surfdir      : path to the directory with surface meshes (use mris_extend for creating intermediate surfaces between white and pial)\n"
					"  --surf_name surfname    : name of a surface file(s) (use * and ? to pass multiple names, maximum %d)\n"
					"                            multiple surface names have to be sorted from wm to pial\n"
					"  --overlay_dir overdir   : path to the directory with surface mesh overlays (use mris_vol2surf to map data onto the surface meshes and create overlays)\n"
					"  --overlay_name overname : name of an overlay file(s) (use * and ? to pass multiple names, maximum %d)\n"
					"                            multiple overlays names have to be sorted from wm to pial\n"
  				"                            corresponding surfaces and overlays must have the same numbers of vertices\n"
  				"  --output_dir outdir     : path to the output directory [default: overdir]\n"
  				"  --output_name outname   : name of a output overlay file [default: based on the name of 1st overlay file]\n"
					"\n"
  				"  --tan-size tansize      : tangential extent of the smoothing kernel [default: 0 = no smoothing, max = %d]\n"
  				"                            depending on the weighting scheme (see --tan-weights):\n"
  				"                            gauss: tansize = FWHM of the 2D gaussian kernel applied in the tangential direction\n"
  				"                            distance: tansize = radius of the neighborhood around the central vertex of the smoothing kernel in tangential direction (number of vertices) \n"
  				"\n"
					"  --rad-size radsize      : radial extent of the intracortical smoothing kernel (number of adjacent meshes) [default: 1 = no smoothing, max = number of input surfaces/overlays]\n"
					"  --rad-start radstart    : starting surface mesh of the intracortical smoothing kernel in the radial direction [default: 0 = white, max = number of input surfaces/overlays--radsize]\n"
  				"\n"
  				"  --tan-weights type      : weighting function for tangential smoothing [default: gauss]\n"
  				"                            gauss = gaussian with FWHM = tansize\n"
  				"                            distance = 1/tansize\n"
  				"  --rad-weights type      : weighting function for radial extent of the kernel: not yet implemented\n"
  				"\n"
  				"  --help                  : prints this info\n"
					"\n"
					"Tangential-only smoothing: \n"
					"  - requires exactly 1 input surface mesh with Nv vertices and exactly 1 corresponding overlay with Nv values\n"
					"  - tansize has to be in range [1, %d]\n"
					"  - one output overlay \n"
					"\n"
					"Intracortical smoothing (radial + tangential): \n"
					"  - requires 2 or more input surfaces spanning from white (=0) to pial (=n)"
					"  - requires 2 or more input overlays - one for each input surface with a corresponding number of vertices\n"
					"  - radsize has to be larger than 2 but smaller than or equal to the number of input surfaces/overlays\n"
					"  - radstart is the index of the starting surface of the intracortical smoothing kernel, where white (=0) & pial (= total number of surfaces-1)\n"
  			  "  - intracortical smoothing kernel will include radsize of adjacent surfaces starting from radstart surface and moving towards pial surface\n"
					"  - tansize has to be in range [0, %d]\n"
  				"\n"
					"If you use mris_smooth_intracortical in your research, please cite:\n"
					"Blazejewska AI, Fischl B, Wald LL, Polimeni JR, Intracortical smoothing of small-voxel fMRI data can provide increased detection power without spatial "
					"resolution losses compared to conventional large-voxel fMRI data.\nNeuroImage 2019. 189:601-614. DOI: 10.1016/j.neuroimage.2019.01.054\n"
					"\n", MAX_SURF, MAX_SURF, MAX_NB, MAX_NB, MAX_NB);
	print_version();
	printf("\n");
	exit(1);
}


/*!
\fn static void print_version(void)
\brief Prints version and exits
*/
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}

/*!
\fn static void print_help(void)
\brief Prints help and exits
*/
static void print_help(void) {
	printf("Smooths data overlaid onto to the cortical surface meshes using kernels for which tangential and radial extent can be specified. \n\n");
	print_usage() ;
	printf("\nWARNING: this program has not yet been tested!\n");
	exit(1) ;
}
