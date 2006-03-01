/**
 * @file   mri_path2label.c
 * @author Kevin Teich
 * @date   $Date: 2006/03/01 23:14:31 $
 * 
 * @brief  Converts scuba's path file format to a label file.
 *
 */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include "version.h"
#include "label.h"
#include "path.h"
#include "error.h"

#ifdef Darwin
#include "getline.h"
#endif

char *Progname = NULL;

static void print_usage (void) ;
static void print_help (void) ;
static int  guess_file_type (char* fname, int* is_path, int *is_label);
static int  convert_path_to_label (char* fname, char* ofname);

struct option long_options[] = 
  {
    {"path2label", 0, 0, 0},
    {"label2path", 0, 0, 0},
    {"help", 0, 0, 0},
    {0, 0, 0, 0}
  };

int main(int argc, char *argv[]) {

  int   nargs;
  int   c;
  int   option_index;
  char* source_file          = NULL;
  char* dest_file            = NULL;
  int   path_to_label  = 0;
  int   label_to_path  = 0;
  int   source_is_label      = 0;
  int   source_is_path       = 0;
  int   err                  = 0;
  FILE* fp                   = NULL;

  nargs = handle_version_option (argc, argv, "$Id: mri_path2label.c,v 1.3 2006/03/01 23:14:31 kteich Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;

  while (1)
    {
      c = getopt_long (argc, argv, "h", long_options, &option_index);
      if (c == -1)
	break;

      switch (c)
	{
	  /* Parsing a long option. */
	case 0:
	  switch (option_index)
	    {
	    case 0: // path2label
	      if (label_to_path)
		{
		  printf ("ERROR: Only specify one option.\n\n");
		  print_usage();
		  exit(0);
		}
	      path_to_label = 1;
	      break;
	    case 1: // label2path
	      if (path_to_label)
		{
		  printf ("ERROR: Only specify one option.\n\n");
		  print_usage();
		  exit(0);
		}
	      label_to_path = 1;
	      break;
	    case 2:
	      print_usage();
	      print_help();
	      exit(0);
	    }
	  break;
	  
	case 'h':
	  print_usage();
	  print_help();
	  exit (0);
	  break;
	}
    }
  
  if ((argc - optind) != 2)
    {
      printf ("ERROR: Please specify an input and output file.\n\n");
      print_usage ();
      exit (0);
    }

  source_file = strdup (argv[optind++]);
  dest_file = strdup (argv[optind++]);

  if (!label_to_path && !path_to_label)
    {
      err = guess_file_type (source_file, &source_is_path, &source_is_label);
      if (err)
	{
	  printf ("ERROR: Couldn't determing source file type.\n");
	  exit (1);
	}

      if (source_is_path)
	path_to_label = 1;
      if (source_is_label)
	label_to_path = 1;
    }

  fp = fopen (dest_file, "w");
  if (NULL == fp)
    {
      printf ("ERROR: Couldn't open dest file for writing.\n");
      exit (1);
    }
  fclose (fp);

  printf ("Converting %s\n", source_file);
  printf ("        to %s\n", dest_file);
  if (path_to_label)
    printf ("Path to label\n");
  else 
    printf ("Label to path\n");
  printf ("\n");

  if (path_to_label)
    convert_path_to_label (source_file, dest_file);

  return 0;
}


static void print_usage(void)
{
  printf ("USAGE: %s [options] input output\n", Progname);
  printf ("\n");
  printf ("   --path2label : will treat input as a path and output a label\n");
  printf ("   --label2path : will treat input as a label and output a path\n");
  printf ("\n");
}


static void print_help(void)
{
  printf (
"  Purpose: Converts a path file to a label or a labe file to a path\n"
"  file. Attempts to guess the correct format by looking for\n"
"  .path and .label suffixes and by looking at the first line in\n"
"  the file. Use one of the options to override this behavior.\n"
"  Will return an error if it cannot guess a format and none is\n"
"  explicitly supplied.\n"
);
  printf ("\n");
}

static int guess_file_type (char* fname, int* is_path, int *is_label)
{
  FILE*  fp         = NULL;
  char*  line       = NULL;
  size_t size       = 1024;
  char*  needle     = NULL;
  int    found      = 0;

  fp = fopen (fname, "r");
  if (NULL == fp)
    {
      printf ("ERROR: Couldn't open %s\n", fname);
      return 1;
    }

  *is_path = 0;
  *is_label = 0;
  found = 0;

  line = (char*) malloc (size);
  
  while (!feof(fp) && !found)
    {
      getline (&line, &size, fp);
      if (line[0] == '#')
	{
	  needle = strstr( line, "label" );
	  if( NULL != needle ) 
	    {
	      *is_label = 1;
	      found = 1;
	      break;
	    }
	  needle = strstr( line, "Path" );
	  if( NULL != needle ) 
	    {
	      *is_path = 0;
	      found = 1;
	      break;
	    }
	}
    }

  fclose( fp );

  free (line);

  if (!found)
    return 1;
  else
    return (ERROR_NONE);
}

static int convert_path_to_label (char* fname, char* ofname)
{
  int num_paths;
  PATH** paths;
  int path_index;
  int pno;
  LABEL* label = NULL;
  int label_vno;

  PathReadMany (fname, &num_paths, &paths);

  label_vno = 0;
  for (path_index = 0; path_index < num_paths; path_index++)
    {
      if (NULL == label)
	{
	  label = LabelAlloc (paths[path_index]->n_points, NULL, NULL);
	}
      else
	{
	  LabelRealloc (label, 
			label->n_points + paths[path_index]->n_points + 1);
	}

      for (pno = 0; pno < paths[path_index]->n_points; pno++)
	{
	  label->lv[label_vno].x = paths[path_index]->points[pno].x;
	  label->lv[label_vno].y = paths[path_index]->points[pno].y;
	  label->lv[label_vno].z = paths[path_index]->points[pno].z;
	  label->lv[label_vno].vno = -1;
	  label_vno++;
	}

      /* Sentinel value. */
      label->lv[label_vno].x = -99999;
      label->lv[label_vno].y = -99999;
      label->lv[label_vno].z = -99999;
      label->lv[label_vno].vno = -1;
      label_vno++;

      PathFree (&paths[path_index]);
    }

  free (paths);

  LabelWrite (label, ofname);

  LabelFree (&label);

  return (ERROR_NONE);
}
