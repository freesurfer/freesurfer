function runlist = fast_runlistfile(runlistfile)
% runlist = fast_runlistfile(runlistfile)

runlist = [];

if(nargin ~= 1)
  msg = 'runlist = fast_runlistfile(runlistfile)';
  qoe(msg);error(msg);
end

fid = fopen(runlistfile);
if(fid == -1)
  msg = sprintf('Could not open %s',runlistfile);
  qoe(msg);error(msg);
end

runid = deblank(fscanf(fid,'%s',1));
while( ~isempty(runid) )
  runlist = strvcat(runlist,runid);
  runid = fscanf(fid,'%s',1);
end

if(isempty(runlist))
  msg = sprintf('ERROR: no runs found in %s\n',runlistfile);
  qoe(msg);error(msg);
end

return;
