function edge = mesh_adjacency(FV)

% modified from mesh_edges, to compute sparse adjacency matrix
% see also triangulation2adjacency from G. Peyré in toolbox_graph
% mesh_edges - Calculate edge lengths of triangulation
% 
% edge = mesh_edges(FV)
% 
% FV.vertices   - vertices of mesh, Nx3 Cartesian XYZ
% FV.faces      - triangulation of vertices
% 
% edge          - edge lengths, indexed by vertex 
%                 number (sparse NxN matrix)
% 
% 
% Licence:  GNU GPL, no implied or express warranties
% History:  07/2002, Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
fprintf('...searching for mesh edges...');

nvertex = size(FV.vertices,1);
nface   = size(FV.faces,1);

% the 'edge' matrix is the connectivity of all vertices
edge = sparse(nvertex,nvertex);

for f = 1:nface,
    
    % compute the length of all triangle edges (Diff is [3x3])
%     Diff = [FV.vertices(FV.faces(f,[1 2 3]),:) - FV.vertices(FV.faces(f,[2 3 1]),:)];
%     Norm = sqrt( sum(Diff.^2, 2) );
    
    edge(FV.faces(f,1),FV.faces(f,2)) = 1;
    edge(FV.faces(f,2),FV.faces(f,3)) = 1;
    edge(FV.faces(f,3),FV.faces(f,1)) = 1;
    
    % make sure that all edges are symmetric
    edge(FV.faces(f,2),FV.faces(f,1)) = 1;
    edge(FV.faces(f,3),FV.faces(f,2)) = 1;
    edge(FV.faces(f,1),FV.faces(f,3)) = 1;
end

t=toc;
fprintf('done (%5.2f sec).\n',t);

return
