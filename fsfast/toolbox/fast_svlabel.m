function err = fast_svlabel(vtxno, xyz, vstat, hdr, filename)
% err = fast_svlabel(vtxno, xyz, vstat, hdr, filename)
%

err = 0;

if(nargin ~= 5)
  msg = 'USAGE: err = fast_svlabel(vtxno, xyz, vstat, hdr, filename)';
  err = 1;
  qoe(msg);error(msg);
end

nvtxs = size(xyz,1);
if(isempty(vtxno)) vtxno = zeros(nvtxs,1); end
if(isempty(vstat)) vstat = zeros(nvtxs,1); end

if(nvtxs  ~= length(vtxno))
  msg = sprintf('nvtxs = %d, != length(vtxno) = %d\n',nvtxs,length(vtxno));
  err = 1;
  qoe(msg);error(msg);
end

if(nvtxs  ~= length(vstat))
  msg = sprintf('nvtxs = %d, != length(vstat) = %d\n',nvtxs,length(vstat));
  err = 1;
  qoe(msg);error(msg);
end

tmp = [vtxno xyz vstat];

fid = fopen(filename,'w');
if(fid == -1)
  msg = sprintf('Could not open %s\n',filename);
  err = 1;
  qoe(msg);error(msg);
end

fprintf(fid,'%s\n',hdr);
fprintf(fid,'%d\n',nvtxs);
fprintf(fid,'%4d   %8.3f  %8.3f  %8.3f     %f\n',tmp'); %'

fclose(fid);

return;
