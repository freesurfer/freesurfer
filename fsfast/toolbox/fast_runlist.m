function runlist = fast_runlist(dirname,runlistfile)
% runlist = fast_runlist(dirname,<runlistfile>)
% runlistfile, if present, should be relative to dirname

runlist = [];

if(nargin ~= 1 & nargin ~= 2)
  msg = 'runlist = fast_runlist(dirname,<runlistfile>)';
  qoe(msg);error(msg);
end

d = dir(dirname);
if(isempty(d))
  msg = sprintf('No runs found in %s\n',dirname);
  fprintf('%s',msg);
  return;
  %qoe(msg); error(msg);
end

if(~exist('runlistfile')) runlistfile = ''; end
if(~isempty(runlistfile))
  rlf = sprintf('%s/%s',dirname,runlistfile);
  runlist = fast_runlistfile(rlf);
  % Should check that they exist too
  return;
end

for n = 1:length(d);
  dname = d(n).name;
  if(length(dname)==3)
    if(~isempty(str2num(dname))) 
      runlist = strvcat(runlist,dname);
    end
  end
end

if(isempty(runlist))
  msg = sprintf('ERROR: no runs found in %s\n',dirname);
  qoe(msg);error(msg);
end

return;
