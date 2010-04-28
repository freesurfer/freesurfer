function err = fast_svpar4(parfile,par4,labels)
% err = fast_svpar4(parfile,par4,labels)

err = 1;
if(nargin ~= 3)
  fprintf('err = fast_svpar4(parfile,par4,labels)\n');
  return;
end

fp = fopen(parfile,'w');
if(fp == -1)
  fprintf('ERROR: opening %s\n',parfile);
  return;
end

nev = size(par4,1);
for nthev = 1:nev
  label = labels(par4(nthev,2)+1,:);
  fprintf(fp,'%9.4f   %2d   %7.3f   %f   %s\n',par4(nthev,1),par4(nthev,2),...
	  par4(nthev,3),par4(nthev,4),label);
end
fclose(fp);

err = 1;
return;
