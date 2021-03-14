

%
% read_csf_patch.m
%
% Original Author: Bruce Fischl
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%


% Should run 'left or right' scripts to set up fname_* variable,s
% and set nfig so that plotting starts in the right location.
fclose('all');
fid=fopen(fname_patch,'rt');
s = fgetl(fid);
s = fgetl(fid);
[poly,count] = sscanf(s,'%d');
numvert = poly(1);
numquad = poly(2);
vertex_coordinates = zeros(3,1);
face_index = zeros(4,1);
vertx_list = zeros(numvert,3);
face_list = zeros(numquad,5);
tic;
for vert = 1:1:(numvert),
    s = fgetl(fid);
    vertx = sscanf(s,'%d');
    s = fgetl(fid);
    vertx_coordinates = sscanf(s,'%f');
    vertx_list(vert,:) = [vertx_coordinates(1:2)' vertx];
end;
toc
tic;
for face = 1:1:(numquad),
    s = fgetl(fid);
    facenum = sscanf(s,'%d');
    s = fgetl(fid);
    face_vertx = sscanf(s,'%d');
    face_list(face,:) = [face_vertx' facenum];
end;
toc
fclose(fid);
full_vertx=zeros(max(max(vertx_list))+1,3);
full_vertx(vertx_list(:,3)+1,1)=vertx_list(:,1);
full_vertx(vertx_list(:,3)+1,2)=vertx_list(:,2);
mesh_real=File2Var(fname_real);
vertx_values=zeros(max(max(vertx_list))+1,1);
vertx_values(vertx_list(:,3)+1)=mesh_real(:,5);

strpwd=pwd;
subplot(2,1,nfig+1);
title([strpwd((length(strpwd)-25):length(strpwd)) ' ' fname_real]);
p_eccen_handle=patch('Vertices',full_vertx,'Faces',face_list(:,1:4)+1,'FaceVertexCData',vertx_values,'FaceColor','interp','EdgeColor','none');
%colormap(rgb(256));
colorbar;

new_val = avg_vertx(vertx_list(:,1),vertx_list(:,2),mesh_real(:,5),4);
new_idx = insert_vals(new_val,vertx_list(:,3));

subplot(2,1,nfig+2);
p_eccen_handle=patch('Vertices',full_vertx,'Faces',face_list(:,1:4)+1,'FaceVertexCData',new_idx,'FaceColor','interp','EdgeColor','none');
title(['Smoothed: 'strpwd((length(strpwd)-25):length(strpwd)) ' ' fname_real]);

