function err = fast_svbslice(y,stem,sliceno,bext,bhdrstr)
% err = fast_svbslice(y,stem,sliceno,<bext>)
%
% size(y) = [rows cols frames]
% sliceno is zero-based
% if bext is not specificed or is null, it is set to bfloat
% bhdrstr is printed to stem.bhdr
%
% if sliceno is < 0 then size(y) = [slices rows cols frames]
% and each slice is saved.
%
% $Id: fast_svbslice.m,v 1.3 2003/07/31 15:12:06 greve Exp $

err = 1;

if(nargin < 2 & nargin > 4)
  fprintf('err = fast_svbslice(y,stem,sliceno,<bext>,<bhdrstr>)\n');
  return;
end

if(exist('bext') ~= 1) bext = ''; end
if(isempty(bext)) bext = 'bfloat'; end

if(exist('bhdrstr') ~= 1) bhdrstr = ''; end

if(strcmp(bext,'bfloat') == 0 & strcmp(bext,'bshort') == 0)
  fprintf('ERROR: bext = %s, must be bfloat or bshort\n',bext);
  return;
end
  
if(sliceno >= 0)
  fname = sprintf('%s_%03d.%s',stem,sliceno,bext);
  fmri_svbfile(y,fname);
else
  nslices = size(y,1);
  for slice = 0:nslices-1
    fname = sprintf('%s_%03d.%s',stem,slice,bext);
    fmri_svbfile(squeeze(y(slice+1,:,:,:)),fname);
  end
end

if(~isempty(bhdrstr))
  % Should make sure nframes is correct %
  bhdrfile = sprintf('%s.bhdr',stem);
  fid = fopen(bhdrfile,'w');
  fprintf(fid,'%s',bhdrstr);
  fclose(fid);
end

err = 0;

return;







