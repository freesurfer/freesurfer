function rtn = dumphdrest(hstem,s,r,c,outfile)
% rtn = dumphdrest(hstem,s,r,c,outfile)
% Example: dumphdrest('h',10,36,48,'h-10-36-48.dat');


%
% dumphdrest.m
%
% Original Author: Doug Greve
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

rtn = 1;

if(nargin ~= 5)
  fprintf('USAGE: r = dumphdrest(hstem,s,r,c,outfile)\n');
  return;
end  

fprintf('Loading data ');tic;
%[havg0 evar0 sxadat] = fast_ldsxabvolume('h');
%[havg0 evar0 sxadat] = fast_ldsxabfile('h_010.bfloat');
infile = sprintf('%s_%03d.bfloat',hstem,s);
[havg0 evar0 sxadat] = fast_ldsxabfile(infile);
fprintf('done %g\n',toc);

%[ns nr nc nh] = size(havg0);
[nr nc nh] = size(havg0);
ns = 1;
nv = ns*nr*nc;

havg = reshape(havg0,[nv nh])'; %'
evar = reshape(evar0,[nv 1])'; %'

dof = sxadat.DOF;
hCM = sxadat.hCovMtx;
dhCM = diag(hCM);

hstderr = sqrt(dhCM * evar);

havg2 = reshape(havg, [sxadat.Nh sxadat.Nnnc ns nr nc]);
hstderr2 = reshape(hstderr, [sxadat.Nh sxadat.Nnnc ns nr nc]);

nn = [1:sxadat.Nh]'; %'
t = sxadat.TER*(nn-1) - sxadat.TPreStim;

s =  1; % Reset slice number because only loaded one slice
%r = 36; 
%c = 48;

%plot(t,0,'b.',t,havg2(:,1,s,r,c),'g',t,havg2(:,2,s,r,c),'r',t,havg2(:,3,s,r,c),'c');
%plot(t,zeros(size(t)),'k-.');
%hold on;
%errorbar(t,havg2(:,1,s,r,c),hstderr2(:,1,s,r,c),'g');
%hold off;

tmp1 = squeeze(havg2(:,:,s,r,c));
tmp2 = squeeze(hstderr2(:,:,s,r,c));

fmt = ['%7.4f '];
fmt = [fmt repmat('%g ',[1 sxadat.Nnnc*2]) ];
fmt = [fmt '\n'];

%outstem = 'tmp';
%outfile = sprintf('%s-s%02d-r%03d-c%03d.dat',outstem,s,r,c);
outfileid = fopen(outfile,'w');
if(outfileid == -1)
  fprintf('ERROR: could not open %s\n',outfile);
  return;
end

fprintf(outfileid,fmt,[t tmp1 tmp2]'); %'
fclose(outfileid);

rtn = 0;

return;
