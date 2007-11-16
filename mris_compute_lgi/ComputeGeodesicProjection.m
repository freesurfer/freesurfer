function [wholepath] = ComputeGeodesicProjection(mesh_total, A, start_vertex, verticeslist)

% --- Step 1: iterative reorganization of the vertices' list 
% (each vertex must appear in the right order so that they can be
% subsequently linked together one by one). 
% This uses mesh_vertex_nearest (D. Weber) and is based on the assumption 
% that verticeslist are relatively regularly situated on a circle. 

% NOTE: sometimes this function can return a small loop instead of the full
% perimeter and it's mandatory to check that with e.g. the size (number of 
% vertices of the perimeter. Just starting the reorglist
% with another start_vertex usually overcomes that problem. 

reorglist=verticeslist(start_vertex);
remaininglist=verticeslist;
remaininglist(start_vertex)=[];

% search the nearest vertex for the first vertex of the list (start_vertex)
[nextindex,nextvalue]=mesh_vertex_nearest(mesh_total.vertices(remaininglist,:),mesh_total.vertices(verticeslist(start_vertex),:));
clear nextvalue; 

% create the reorglist (reorganized list)
reorglist=[reorglist remaininglist(nextindex)];

% delete this vertex from the list of vertices (i.e. from the pool of
% vertices to be reorganized)
remaininglist(nextindex)=[];

% repeat that one more time so that we are sufficiently far from the 
% start_vertex to add it to the pool of vertices to be reorganized. 
% Then continue iteratively until we close the loop by finding the start_vertex. 
[nextindex,nextvalue]=mesh_vertex_nearest(mesh_total.vertices(remaininglist,:),mesh_total.vertices(reorglist(2),:));
clear nextvalue;
reorglist=[reorglist remaininglist(nextindex)];
remaininglist(nextindex)=[];

% add the start_vertex
remaininglist= [remaininglist verticeslist(start_vertex)];

% continue to reorganize until start_vertex
for z= 3: size(verticeslist,2)     
    [nextindex,nextvalue]=mesh_vertex_nearest(mesh_total.vertices(remaininglist,:),mesh_total.vertices(reorglist(z),:));
    clear nextvalue;
    reorglist=[reorglist remaininglist(nextindex)];
    if remaininglist(nextindex) == verticeslist(start_vertex), break, end
    remaininglist(nextindex)=[];
end

clear remaininglist;
clear nextindex;
clear nextvalue;


% ---- Step 2: find the shortest geodesic path to link each vertex of the 
% reorglist. This uses dijk, a function by Michael G. Kay (Matlog toolbox);
% which itself requires the adjacency matrix (A).
% The whole perimeter of the pial ROI is stored in the vector "wholepath".

k=size(reorglist,2)-1;
wholepath=reorglist(1);

for q=1:k
    [D,P]=dijk(A,reorglist(q),reorglist(q+1));
    wholepath=[wholepath P(2:end)];
end

