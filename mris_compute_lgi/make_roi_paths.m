function make_roi_paths(pial, outersmoothed, radius, stepsize, outdir, flagf)

% Original Author: Marie Schaer
% Date: 2007/11/15
%
% This function creates multiple circular regions of interest on the outer
% surface, measures their area (saved in matlab format) and save path file
% containing a set of
% corresponding points on the pial surface. The path files are required as
% an input to the mri_path2label function
%
% Parameters:
% radius: radius of the circular region of interest on the outer surface.
% Most appropriate values between 20 and 25 mm, further documentation
% available in the validation publication: "A surface-based approach to
% quantify local cortical gyrification" by Schaer M, Bach Cuadra M,
% Tamarit L, Eliez E, Thiran JP, IEEE Transactions on Medical Imaging, 2007
%
% Utilities to transfer points on the pial surface calls "mesh_vertex_nearest",
% a function written by Darren Weber. Utilities to compute matrix adjacency
% for the pial surface calls "mesh_adjacency", a program first written by
% Darren Weber, then modified by G. Peyre (source: toolbox_graph on
% mathworks)
%
% Example: make_roi_paths('lh.pial','lh.outer-pial-smoothed',25,100,/tmp)

t0 = cputime;

disp('loading datas ...')
[mesh_total.vertices, mesh_total.faces]=freesurfer_read_surf(pial);
[mesh_outer.vertices, mesh_outer.faces]=freesurfer_read_surf(outersmoothed);

disp('preparing outer mesh structure ...')
mesh_outer = createMeshFacesOfVertex (mesh_outer.vertices, mesh_outer.faces)

disp('preparing pial mesh structure ...')
A = mesh_adjacency(mesh_total);

% Keep track of the area for each individual region of interest on the
% outer surface:
nbVertices = size(mesh_outer.vertices,1);
areasOuterROI = zeros(nbVertices,1);
areasOuterROI = sparse (areasOuterROI);
clear nbVertices


% Circular regions of interest are defined each 100 vertex on the outer
% surface. Due to the high resolution of the outer mesh (average face area
% of ~0.3mm2), calculations are computed each 1 on 100 vertex, to avoid
% high redundancies and optimize calculation time.

for iV = 1:stepsize: length(mesh_outer.vertices)

    disp(['... creating path file for vertex ',num2str(iV),' / ',num2str(length(mesh_outer.vertices))])

    % ------ Part 1: find the ROI on the enveloppe --------------

    [verticesInROI, facesInROI] = getVerticesAndFacesInSphere(mesh_outer, iV, radius);

    % In rare case, the intersection of the enveloppe with the sphere result
    % in two regions of interest: the one that we expect, circular and
    % centered in iV; and an aberrant one at the bottom of the sphere
    % (e.g. near the superior sagittal vault). The next step is used to
    % control for that and keep only the radial region centered at iV.
    facesListOuterGeo = MakeGeodesicOuterROI (mesh_outer, iV, facesInROI);

    % Find the perimeter of the outer ROI
    verticesOuterGeo=unique(facesListOuterGeo(:));
    perim=setdiff(verticesOuterGeo,verticesInROI);

    % Measure its area
    areaOuterROI = getFacesArea(mesh_outer,facesListOuterGeo);

    areasOuterROI(iV) = areaOuterROI;

    if (mod(iV,50*stepsize) == 1)
        fname = [outdir '/' pial '.area_intermed.mat'];
        save(fname, 'areasOuterROI') ;
        disp(['area file for outer ROIs saved at ',num2str(iV)])
    end

    clear verticesInROI facesInROI facesListOuterGeo verticesOuterGeo areaOuterROI


    % ----- Part2: Define the corresponding set of points on the pial surface:
    % Transfer a few points of the perimeter of the hulls ROI from
    % the envelope to the pial surface, and save them as path file.
    step = 7;
    [verticeslist] = SearchProjectionOnPial(mesh_total, mesh_outer, perim, step);

    % reorganize the set of points in the right order to input them to
    % mri_path2label
    reorglist = reorganize_verticeslist (mesh_total, A, mesh_outer, perim, verticeslist, step);
    lindex = reorglist';
    lxyz = mesh_total.vertices(lindex,:); lindex= lindex-1;


    p = sprintf ('%06d', iV);
    write_path(lxyz, lindex,[outdir '/' pial '.' p '.path']);

    clear verticeslist lindex lxyz

end

disp(['saving area file in matlab format... '])
fname = [outdir '/' pial '.outer_ROIs_area.mat'];
save(fname, 'areasOuterROI')

disp('DONE.')
deltaT = cputime - t0

% indicate successful completion by deleting flagfile
delete(flagf);
