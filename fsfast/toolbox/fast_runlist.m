function runlist = fast_runlist(dirname)
% runlist = fast_runlist(dirname)

runlist = [];

if(nargin ~= 1)
  msg = 'runlist = fast_runlist(dirname)';
  qoe(msg);error(msg);
end

d = dir(dirname);
if(isempty(d))
  msg = sprintf('No runs found in %s\n',dirname);
  qoe(msg); error(msg);
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
