#ifndef PATHS_INCLUDED
#define PATHS_INCLUDED

typedef struct
{
  float x;
  float y;
  float z;
} PATH_POINT;

typedef struct
{
  int n_points;         /* number of points in this path. */
  char name[100];       /* original file name */
  PATH_POINT *points;   /* array of size n_points */
} PATH;


int    PathReadMany (char *fname, int *num_read, PATH ***paths);
int    PathWriteMany (char *fname, int num_paths, PATH **paths);
PATH*  PathAlloc (int n_points, char* name);

#endif
