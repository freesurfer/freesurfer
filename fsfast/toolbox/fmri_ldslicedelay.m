function [slc, slcdel] = fmri_ldslicedelay(fname)
% [slc slcdel] = 

if(nargin ~= 1)
  msg = 'USAGE: [slc slcdel] = fmri_ldslicedelay(fname)';
  qoe(msg);error(msg);
end

fname = deblank(fname);

fid = fopen(fname,'r');
if(fid == -1)
  msg = sprintf('Could not open %s',fname);
  qoe(msg);error(msg);
end

tmp = fscanf(fid,'%f');
fclose(fid);

nslices = length(tmp)/2;
tmp = reshape(tmp, [2 nslices])'; %'
slc    = tmp(:,1);
slcdel = tmp(:,2);

%% resort by slice number %%
[y i] = sort(slc);
slc    = y;
slcdel = slcdel(i);

return;
