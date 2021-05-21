function spm_to_bshort(iname)
% This function reads an SPM image volume xxx.img and writes it out 
% as a series of bshort slices xxx_num.img xxx_num.hdr.
% FORMAT spm_to_bshort(iname);
% iname   - input image name
%______________________________________________________________________
% Writes out transaxially sliced SPM image volume into a series of
% slices *.bshort. 
% (The bshort format is for MGH flattening compatibility)
%___________________________________________________________________________
% @(#)spm_to_bshort.m            99/31/03                 Chloe Hutton


%
% spm_to_bshort.m
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

%
% Use user interface if no input
%-----------------------------------------------------------------------
if nargin==0
   a=spm_get(1,'.img','select');
else
  a=iname;
end;  
[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(a);

Vm=spm_map(a,[DIM VOX 1 TYPE OFFSET]);
D=spm_matrix([0 0 0]);
fprintf('Writing slices ...');
TYPE=4;

% Write out slices
%--------------------------------------------------------------------------
for i=1:DIM(3)
   if i<=10
      num=sprintf('_00%d',i-1);
   elseif i<=100
      num=sprintf('_0%d',i-1);
   elseif i<=1000
      num=sprintf('_%d',i-1);
   end

   dfname=deblank(spm_str_manip(a,'r'));
   name=[dfname,num,'.bshort'];
   hname=[dfname,num,'.hdr'];

   fid=fopen(name,'w');
   D(3,4)=-i;
   slice=spm_slice_vol(Vm,inv(D),[DIM(1) DIM(2)],1);
   fwrite(fid,slice,spm_type(TYPE));
   fclose(fid);
   fid=fopen(hname,'w');
   fprintf(fid,'%d %d 1 0\n\n',DIM(1),DIM(2));
   fclose(fid);
end






