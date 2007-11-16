function [nearestIndex,nearestValues] = mesh_vertex_nearest(vertices,points)

% mesh_vertex_nearest - find nearest vertices to specified points
%
% Usage: [nearestIndex,nearestValues] = mesh_vertex_nearest(vertices,points)
%
% vertices is a Vx3 matrix of 3D Cartesian coordinates.
% points is a Px3 matrix of 3D Cartesian coordinates.  These points need not
% be among the vertices, but they are somewhere near to particular points 
% in the vertices cloud.  The function finds just one of the nearest 
% vertices in the cloud for each of these points.
%
% nearestIndex is the indices into vertices nearest to points
% nearestValues is the coordinates for nearestIndex
%
% This function is just a wrapper for dsearchn.
%

% $Revision: 1.1 $ $Date: 2007/11/16 02:54:47 $

% Author:   Darren.Weber_at_radiology.ucsf.edu
% Created:  23 May 2004
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fprintf('...Finding nearest vertices to points.\n');


nearestIndex = dsearchn(vertices,points);
nearestValues = vertices(nearestIndex,:);

% for p = 1:length(points),
%     tmp = repmat(points(p,:),length(vertices),1);
%     [minValue, nearestIndex] = min(abs(vertices - tmp));
% end
% nearVertex = vertices(nearestIndex,:);





return
