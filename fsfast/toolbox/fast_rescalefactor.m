function [RescaleFactor, MeanVal] = fast_rescalefactor(MeanValFile, RescaleTarget)
% fast_rescalefactor(MeanValFile, RescaleTarget)

if(nargin ~= 2)
  msg = 'USAGE: [factor meanval] = fast_rescalefactor(MeanValFile, RescaleTarget)';
  qoe(msg);error(msg);
end

fid = fopen(MeanValFile);
if(fid == -1)
  msg = sprintf('ERROR: could not open %s',MeanValFile);
  qoe(msg);error(msg);
end

MeanVal = fscanf(fid,'%f');
RescaleFactor = RescaleTarget/MeanVal;

fclose(fid);

return;






