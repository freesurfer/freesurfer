function r = fast_fileexists(filename)
% r = fast_fileexists(filename)
% 1 if it exists and is readable , 0 if not

if(nargin ~= 1) 
  msg = 'USAGE: r = fast_fileexists(filename)';
  qoe(msg); error(msg);
end

fid = fopen(filename,'r');
if(fid == -1 ) r = 0;
else
  r = 1;
  fclose(fid);
end

return;
