#ifndef GEODESIC_DLL_HEADER_HPP_234232
#define GEODESIC_DLL_HEADER_HPP_234232

#ifndef GEODESIC_DLL_IMPORT 
#define GEODESIC_DLL_IMPORT
#endif

#ifdef __cplusplus
extern "C" {
#endif

long new_mesh(long num_points,		//creates new mesh
								  double* points,	
								  long num_triangles,
								  long* triangles, 
								  long* num_edges, 
								  double** edges);

long new_algorithm(long mesh_id,	//creates a geodesic algorithm for a given mesh
						               long type,
									   long subdivision);

void delete_algorithm(long id);

void delete_mesh(long id);			//delete mesh and all associated algorithms

void propagate(long algorithm_id,		//compute distance field for given source points
									double* source_points,	
									long num_sources,
									double* stop_points,	//limitations on distance field propagation
									long num_stop_points,
									double max_propagation_distance);

long trace_back(long algorithm_id,		//using procomputed distance field, compute a shortest path from destination to the closest source
									double* destination,
									double** path);

long distance_and_source(long algorithm_id,		//quickly find what source this point belongs to and what is the distance to this source
											 double* destination,			
											 double* best_source_distance);

long distance_and_source_for_all_vertices(long algorithm_id,	//same idea as in the previous function
															  double** distances,	//list distance/source info for all vertices of the mesh
															  long** sources);

#ifdef __cplusplus
}
#endif

#endif
