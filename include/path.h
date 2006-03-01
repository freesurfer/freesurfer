#ifndef PATHS_INCLUDED
#define PATHS_INCLUDED

/* One of these for every point in the path. */
typedef struct
{
  float x;
  float y;
  float z;
} PATH_POINT;

/* Main path struct. */
typedef struct
{
  int n_points;         /* number of points in this path. */
  char name[100];       /* original file name */
  PATH_POINT *points;   /* array of size n_points */
} PATH;


/* Read in multiple paths from a path file. On return, num_read will
   point to the number of paths in the paths variable. Pass in the
   address to a PATH** variable. On return it will point to a list of
   allocated PATH objects. */
int    PathReadMany (char *fname, int *num_read, PATH ***paths);

/* Write multiple paths to a file. paths should be an array of
   pointers to PATH structs, each one with a valid PATH
   object. num_paths should be the number of paths in that array. */
int    PathWriteMany (char *fname, int num_paths, PATH **paths);

/* Allocate and free paths. */
PATH*  PathAlloc (int n_points, char* name);
int    PathFree (PATH** path);

#endif
