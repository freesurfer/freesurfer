function spm_to_bfloat(iname)
% This function reads a single byte SPM image and writes it out
% as a series of float slices called actfloat_%d.bfloat
% FORMAT spm_to_bfloat(iname);
% iname   - input image name
%______________________________________________________________________
% Writes out transaxially sliced SPM activation images into a series of
% slices actfloat_%d.img. 
% (The actfloat format is for MGH flattening compatibility)
%___________________________________________________________________________
% @(#)spm_to_bfloat.m            99/31/03                 Chloe Hutton


%
% spm_to_bfloat.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:34 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

% Use user interface if no input
%-----------------------------------------------------------------------
if nargin==0
   a=spm_get(1,'.img','select');
else
   a=iname;
end;  
[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(a);

% Output slices are written to actfloat_%d
%----------------------------------------------------------------------
pathname=deblank(spm_str_manip(a,'h'));
dfname=[pathname,'/','actfloat_'];

Vm=spm_map(a,[DIM VOX SCALE TYPE OFFSET]);
D=spm_matrix([0 0 0]);

for i=1:DIM(3)

   if i<=10
      num=sprintf('00%d',i-1);
   elseif i<=100
      num=sprintf('0%d',i-1);
   elseif i<=1000
      num=sprintf('%d',i-1);
   end
   name=[dfname,num,'.bfloat'];
   hname=[dfname,num,'.hdr'];
   fid=fopen(name,'w');

   % read in each slice and save as a float slice
   %------------------------------------------------------------------
   D(3,4)=-i;
   slice=spm_slice_vol(Vm,inv(D),[DIM(1) DIM(2)],1);
   TYPE=16;
   fwrite(fid,slice,spm_type(TYPE));
   fclose(fid);
   hfid=fopen(hname,'w');
   fprintf(hfid,'%d %d 1 0',DIM(1),DIM(2));
   fclose(hfid);
end






