function find_corresponding_center_FSformat (pial, outersmoothed, outdir)

% Original Author: Marie Schaer
% Date: 2007/11/14
% 
% This function takes as an input the pial mesh and its envelope (the outer
% mesh), and is used to isolate a set of corresponding points of both 
% meshes. In the % lGI computation process, it is used to define the seed 
% point for the obtention of the pial regions of interests (by filling 
% until a closed perimeter).
% 
% The function uses "mesh_vertex_nearest" and "read_surf"
% from Darren Weber's bioelectromagnetism toolbox
%
% Input:
%  ?h.pial surface from FreeSurfer
%  ?h.outer-pial from make_outer_surface
% 
% Output:
%  c_* asci files containing the index of the mesh_pial vertex
%  corresponding to the vertex number * of the outer mesh,
%  and a file named 'center.vertices' which contains selected vertices

t0 = cputime;    

[mesh_pial.vertices, mesh_pial.faces] = read_surf(pial);
[mesh_outer.vertices, mesh_outer.faces] = read_surf(outersmoothed);

fidv = fopen([outdir '/' pial '.center.vertices'], 'w') ;

for iV = 1 :100: length(mesh_outer.vertices) 
        
    centerSeed=mesh_vertex_nearest(mesh_pial.vertices,mesh_outer.vertices(iV,:));
    centerSeed = centerSeed - 1; %FreeSurfer-s vertex are 0-based
    p = sprintf ('%06d', iV);
      
    fid = fopen([outdir '/' pial '.' p '.center'], 'w') ;
    fprintf(fid,'%d\n',centerSeed);
    fclose(fid) ;

    fprintf(fidv,'%s\n',p);
end
    
fclose(fidv) ;

deltaT = cputime - t0
