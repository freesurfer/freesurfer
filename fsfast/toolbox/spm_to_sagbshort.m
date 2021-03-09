function spm_to_sagbshort(iname)
% Writes out a transaxially sliced image volume into a series of sagittal 
% slices with the .bshort suffix.
% FORMAT spm_to_bshort(iname);
% iname   - input image filename
%___________________________________________________________________________
% Writes out transaxially sliced SPM image xxxx.img into series of sagittal 
% slices xxxx_num.bshort each with an associated header file xxxx_num.hdr
% (The bshort format is for MGH flattening compatibility)
%___________________________________________________________________________
% @(#)spm_to_sagbshort.m            99/31/03                 Chloe Hutton


%
% spm_to_sagbshort.m
%
% Original Author: Doug Greve
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

if nargin==0
   VFn=spm_get(1,'.img','select structural image');
else
   VFn=iname;
end;

[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(VFn);
VF=spm_map(VFn,[DIM VOX 1 TYPE OFFSET]);

% this transform converts transverse to coronal
%---------------------------------------------------------------------------
M=[0 1 0 0;0 0 1 0;1 0 0 0;0 0 0 1];
nDIM=[DIM(2) DIM(3) DIM(1)];
nVOX=[VOX(2) VOX(3) VOX(1)];
nORIGIN=[ORIGIN(2) ORIGIN(3) ORIGIN(1)];

% Create temporary file to hold sagittal image volume
%---------------------------------------------------------------------------
iname=deblank(spm_str_manip(VFn,'r'));
dfname=deblank(spm_str_manip(VFn,'h'));
name=[dfname,'/sag_tmp.img'];
hname=[dfname,'/sag_tmp.hdr'];
fprintf('Converting to sagittal ...');
fid=fopen(name,'w');
colormap(gray);
for i=1:nDIM(3)
   M(3,4)=-i;
   slice=spm_slice_vol(VF,inv(M),[nDIM(1) nDIM(2)],-1);
   fwrite(fid,slice,spm_type(TYPE)); 
end
fclose(fid);
spm_hwrite(hname,nDIM,nVOX,SCALE,TYPE,OFFSET,nORIGIN,DESCRIP);
fprintf('\n');

% The second part of this script reslices the sagittal image volume into 
% as series of 256 slices (*.bshort) each with  256*256 bytes.
%--------------------------------------------------------------------------
VFn=name;
[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(VFn);
VF=spm_map(VFn);

% Calculate number of blank spaces so that the image origin is centred 
% in the 256x256 bytes (to keep the Talairach origin in the correct place).
%--------------------------------------------------------------------------
blankx=128-ORIGIN(1);
blanky=128-ORIGIN(2);
blankz=128-ORIGIN(3);

% Set the image type to 4 (bshort)
%------------------------------------------------------------------------- 
TYPE=4;

fprintf('Writing bshorts ...');

% Write out initial slices as blanks if necessary
%------------------------------------------------------------------------- 
for i=1:blankz
   if i<10
   	outfile=sprintf('%s_00%d.bshort',iname,i);
        outhdr=sprintf('%s_00%d.hdr',iname,i);
   elseif i<100
        outfile=sprintf('%s_0%d.bshort',iname,i);
        outhdr=sprintf('%s_0%d.hdr',iname,i);
   else
        outfile=sprintf('%s_%d.bshort',iname,i);
        outhdr=sprintf('%s_%d.hdr',iname,i);
   end   
   fid=fopen(outfile,'w'); 
   nslice=zeros(256,256); 
   slice=nslice;
   fwrite(fid,slice,spm_type(TYPE));
   fclose(fid);
   fid=fopen(outhdr,'w');
   fprintf(fid,'%d %d 1 0\n\n',DIM(1),DIM(2));
   fclose(fid);
end

% Write out slices, padding zeros to make 256*256
%--------------------------------------------------------------------------
D=spm_matrix([0 0 0]);
for i=1:DIM(3)
   ii=i+blankz;
   if ii<10
   	outfile=sprintf('%s_00%d.bshort',iname,ii);
        outhdr=sprintf('%s_00%d.hdr',iname,ii);
   elseif ii<100
        outfile=sprintf('%s_0%d.bshort',iname,ii);
        outhdr=sprintf('%s_0%d.hdr',iname,ii);
   else
        outfile=sprintf('%s_%d.bshort',iname,ii);
        outhdr=sprintf('%s_%d.hdr',iname,ii);
   end   

   fid=fopen(outfile,'w');
   D(3,4)=-i; 
   nslice=zeros(256,256); 
   slice=spm_slice_vol(VF,inv(D),[DIM(1) DIM(2)],1);
   nslice(blankx+1:blankx+DIM(1),blanky+1:blanky+DIM(2))=slice;
   slice=rot90(nslice,2);
   fwrite(fid,slice,spm_type(TYPE));
   fclose(fid);
   fid=fopen(outhdr,'w');
   fprintf(fid,'%d %d 1 0\n\n',DIM(1),DIM(2));
   fclose(fid);
end

% Add any necessary blanks slices at end
%-------------------------------------------------------------------------
for i=DIM(3)+blankz+1:256
   if i<10
   	outfile=sprintf('%s_00%d.bshort',iname,i);
        outhdr=sprintf('%s_00%d.hdr',iname,i);
   elseif i<100
        outfile=sprintf('%s_0%d.bshort',iname,i);
        outhdr=sprintf('%s_0%d.hdr',iname,i);
   else
        outfile=sprintf('%s_%d.bshort',iname,i);
        outhdr=sprintf('%s_%d.hdr',iname,i);
   end   
   fid=fopen(outfile,'w'); 
   nslice=zeros(256,256); 
   slice=nslice;
   fwrite(fid,slice,spm_type(TYPE));
   fclose(fid);
   fid=fopen(outhdr,'w');
   fprintf(fid,'%d %d 1 0\n\n',DIM(1),DIM(2));
   fclose(fid);
end
fprintf('\n');
eval(['!\rm ' name]);
eval(['!\rm ' hname]);















