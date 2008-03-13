//MARS Project (Thomas Yeo and Mert Sabuncu), MIT, CSAIL (c) 2006-2008

#ifndef MARS_DT_Boundary

/*
  sv: int32 array of size N vertices, with the entries corresponding to your boundary marked as 1, and
the remaining entries marked as 0.
numVerts: int32 variable, with value = N
vertNbors: M x N int32 matrix. where the (i, j) entry corresponds to vertex number of the ith neighbor
of vertex j. If vertex j has less than M neighbors, for example 4 neighbors, then entries (i > 4, j)
are set to 0.
maxNeighbors: int32 variable with value = M
costNbors: M x N double matrix, where the (i, j) entry corresponds to distance between ith neighbor of
vertex j and vertex j. Like vertNbors, if the ith neighbor doesn't exist, then the entry is set to 0.
final_cost:  array of size N of type double. MARS_DT_Boundary will save the output into final_cost.
*/

void MARS_DT_Boundary(int *sv, int numVerts, int maxNeighbors, int *vertNbors, double *costNbors, double *final_cost);

/*
  convert from a 2d coord pair to a 1d index (e.g. index = index_2D_array(vno, j, max_nbrs))
*/
int index_2D_array(int row, int col, int num_rows);


#endif
