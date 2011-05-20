function compute_lgi (pial, outersmoothed, fnormal, nsteps, radius, outdir)

% Original Author: Marie Schaer
% Date: 2007/11/15
% 
% This function computes final local GI calculations, based on previously
% measured area for outer surfaces ROIs (outer_ROIs_area.mat) and on their
% corresponding regions of interest on the pial surfaces (stored in 
% ROI_*.label files). 
%
% Input:
%  "pial" is a pial surface in FreeSurfer format
%  "outersmoothed" is the outer smoothed envelope in FreeSurfer format
%  "fnormal" is the ascii file containing info on the normals to the outer
% smoothed surface (as given by mris_convert -n ) 
%  "nsteps" is the interval which was previously used for lGI calculations
%  "radius" is the radius which was previously used for lGI calculations
% 
% Output:
% lGI scalar values sampled on the pial surface and stored in lGI.mgh files
% that can be given as an input to mri_glmfit for further statistical
% analyses on cortical gyrification differences between groups.
% 
% Example: compute_lgi ('lh.pial', 'lh.outer-smoothed', 'lh.outer-smoothed_n.asc', 100 , 25)


t0 = cputime;
    
disp('loading datas ...')
[mesh_total.vertices, mesh_total.faces]=freesurfer_read_surf(pial);
[mesh_outer.vertices, mesh_outer.faces]=freesurfer_read_surf(outersmoothed);

disp('preparing outer mesh structure ...')
mesh_outer = createMeshFacesOfVertex (mesh_outer.vertices, mesh_outer.faces);

disp('reading normals info ...')
vertexNormal = read_normals (mesh_outer, fnormal);
mesh_outer.vertexNormal = vertexNormal;
mesh_outer

disp('preparing pial mesh structure ...')
mesh_total = createMeshFacesOfVertex (mesh_total.vertices, mesh_total.faces)
A = mesh_adjacency(mesh_total);

% vertexRatio will contain lGI measurements for the outer mesh
vertexRatio = zeros(size(mesh_outer.vertices,1),1); 
vertexRatio = sparse (vertexRatio); 
% mt* variables contain info sampled on the pial surface
mtTotalWeightDegressive = zeros(size(mesh_total.vertices,1),1);
mtTotalWeightDegressive = sparse (mtTotalWeightDegressive);
mtTotalRatioDegressive = zeros(size(mesh_total.vertices,1),1);
mtTotalRatioDegressive = sparse (mtTotalRatioDegressive);
    
disp('loading outer area information...')
fname=[outdir '/' pial '.outer_ROIs_area.mat'];
load(fname,'areasOuterROI');


for iV = 1 : nsteps: length(mesh_outer.vertices)
    
    p = sprintf ('%06d', iV);
    
    [ VerticesInPialROI ] = read_ROIlabel ([outdir '/' pial '.' p '.label']); %read_ROIlabel is based on read_label but do not require the subject's name as an input
        
    FacesInPialROI = [];
    for j = VerticesInPialROI'
        FacesInPialROI = [FacesInPialROI, mesh_total.facesOfVertex(j).faceList];
    end
    FacesInPialROI = unique(FacesInPialROI);
    ListFacesInPialROI = mesh_total.faces(FacesInPialROI,:);
    areaPialROI = getFacesArea(mesh_total,ListFacesInPialROI);
    lGI = areaPialROI/areasOuterROI(iV);
   
          
    if lGI > 9 
        disp(['... remeasuring lGI value for vertex iV = ', num2str(iV), '. It may take a few minutes.'])
        areaPialROI = getMeshArea(mesh_total) - areaPialROI ; % see the corresponding label, which is probably inverted
        if areaPialROI < areasOuterROI(iV)
            disp(['WARNING -- Problem for vertex iV = ', num2str(iV), ', lGI value is aberrantly high (lGI=', num2str(lGI), ')...'])
            disp([ '...lGI computation will be stopped. This may be caused by topological defects, check mris_euler_number on the pial surface.'])
            return
        end
        lGI = areaPialROI/areasOuterROI(iV);
        
        if lGI > 9
            disp(['WARNING -- Problem for vertex iV = ', num2str(iV), ', lGI value is aberrantly high (lGI=', num2str(lGI), ')...'])
            disp([ '...lGI computation will be stopped. This may be caused by topological defects, check mris_euler_number on the pial surface.'])
            return
        end
    end

    disp(['lGI for vertex number ',num2str(iV),' of the outer mesh is ',num2str(lGI)])

    vertexRatio(iV) = lGI; % In case the lGI was wrongly computed in both freesurfer and matlab it will crash here...
    clear lGI;
        

        % ----Step 3: Propagate the lGI value back on the pial surface
        VerticesInPialROI=VerticesInPialROI';
        for iVc = VerticesInPialROI
            vAx = transVertexToNormalAxisBase(mesh_outer, iV, mesh_total.vertices(iVc,:));
            distance = norm(vAx(1:2));
            clear vAx;
            weight = 1/(distance+1); % degressive weighting scheme, where the weight is equivalent to 1 if the pial vertex is just on the main axis of the iV normal, and decreases exponentially with the distance to this axis
          
            mtTotalWeightDegressive(iVc) = mtTotalWeightDegressive(iVc) + weight;
            mtTotalRatioDegressive(iVc) = mtTotalRatioDegressive(iVc) + weight*vertexRatio(iV);           
        end

        clear VerticesInPialROI;
        clear weight;
        clear distance;
               
        if (mod(iV,100 * nsteps) == 1)
%            save ratio_intermed.mat mtTotalWeightDegressive mtTotalRatioDegressive vertexRatio ; 
            disp(['ratio saved at ',num2str(iV)])
        end
            
    end
 
%    save ratio_intermed.mat mtTotalWeightDegressive mtTotalRatioDegressive vertexRatio  ;
    
    
fprintf('.. Almost finished .. now back-propagating values on the folded mesh ..')
    
mtActualRatioDegressive = zeros(size(mesh_total.vertices,1),1);
 
% calculate the actual lGI ratios value on the pial mesh
for iV = 1:size(mesh_total.vertices,1)
    if (mtTotalWeightDegressive(iV) > 0)
        mtActualRatioDegressive(iV) = mtTotalRatioDegressive(iV) / mtTotalWeightDegressive(iV);
    else
        mtActualRatioDegressive(iV) = 0;
    end
end

fname=[outdir '/' pial '.outer_lGI.mat'];
save(fname,'vertexRatio');
fname=[outdir '/' pial '.pial_lGI.mat'];
save(fname,'mtActualRatioDegressive');

average_lGI_over_the_hemisphere = mean (mtActualRatioDegressive)

fname=[outdir '/' pial '_lgi.asc'];
write_lgi(mesh_total, mtActualRatioDegressive, fname);

disp('DONE.')
deltaT = cputime - t0
