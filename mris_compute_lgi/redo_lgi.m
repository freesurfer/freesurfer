function [correct_lGI] = redo_lgi(mesh_total, A, mesh_outer, iV, radius, step)

        % recompute previous steps until it get the path file
        [verticesInROI, facesInROI] = getVerticesAndFacesInSphere(mesh_outer, iV, radius);
        facesListOuterGeo = MakeGeodesicOuterROI (mesh_outer, iV, facesInROI);
        verticesOuterGeo=unique(facesListOuterGeo(:));
        perim=setdiff(verticesOuterGeo,verticesInROI);
        areaOuterROI = getFacesArea(mesh_outer,facesListOuterGeo);
        clear verticesInROI facesInROI facesListOuterGeo verticesOuterGeo
        [verticeslist] = SearchProjectionOnPial(mesh_total, mesh_outer, perim, step); % verticeslist correspond to the (not reorganized) set of points which was saved as path_file
        
        % Then, this set of vertices is linked to form a fully continuous
        % perimeter on the pial surface. If two points are situated on
        % different gyri, they are linked through the deepness of the sulci
        % between them. 
        % The function "ComputeGeodesicProjection" calls "mesh_vertex_nearest" (D.Weber)
        % and "dijk" (a function written by Michael G. Kay as a part of the
        % Matlog toolbox).
        pass_number = 1;
        start_vertex = 1;
        disp('starting geodesic region from first vertex...')
        [wholepath] = ComputeGeodesicProjection(mesh_total, A, start_vertex, verticeslist);
        start_vertex=2;
        while ( wholepath(end) ~= wholepath(1) ) | ( size(wholepath,2) < 4*(size(verticeslist,2)) )  
             if start_vertex>=size(verticeslist,2) & (mod(pass_number,2) == 0);
                 step = step +1; % limit redudancies in the path file usually correct the trouble
                 [verticeslist]=SearchProjectionOnPial(mesh_total,m,perim,step);
                 start_vertex=1;
                 pass_number=pass_number+1;
                 disp(['starting geodesic region from first vertex... / step ',num2str(pass_number)])
             end
             fprintf('...bad geodesic region ...recomputing with another startpoint: ')
             disp(num2str(start_vertex))
             [wholepath]=ComputeGeodesicProjection(mesh_total,A,start_vertex,verticeslist);             
             start_vertex=start_vertex+2;   
        end
        disp('...good geodesic region.. propagating...')
        clear start_vertex;
        clear verticeslist;
        
        % At this stage, we have a good perimeter for the pial ROI. The
        % circular ROI is then "filled" with the full geodesical region
        % contained within this perimeter. 
             
        [ListVerticesInPialROI,ListFacesInPialROI]=PropagateGeodesic(mesh_total,mesh_outer,iV,wholepath);
        clear wholepath;
        
        areaMT = getFacesArea(mesh_total,ListFacesInPialROI)
        clear ListFacesInPialROI;
       
        % local Gyrification Index is calculated 
        lGI = areaMT / areaOuterROI;
        
        if lGI <= 7
           correct_lGI = lGI;
        end
            
        