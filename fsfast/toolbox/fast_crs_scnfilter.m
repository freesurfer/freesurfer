function Fcrs = fast_crs_scnfilter(c,r,s,V)
%
% Fcrs = fast_crs_scnfilter(c,r,s,V)
%
% Returns the spatial noise cancellation filter kernal
% for a single point (c,r,s) in the volume. crs are 1-based.
%
% V are the spatial eigenvectors. V is Ns x Nr x Nc x Nev.
%


%
% fast_crs_scnfilter.m
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


Fcrs = [];

if(nargin ~= 4)
  fprintf('USAGE: Fcrs = fast_crs_scnfilter(c,r,s,V)\n');
  return;
end

vdim = size(V);
Ns  = vdim(1);
Nr  = vdim(2);
Nc  = vdim(3);
Nv  = Ns*Nr*Nc;
Nev = vdim(4);

if(s < 1 | s > Ns)
  fprintf('ERROR: s=%d, out of range [1 %d]\n',s,Ns);
  return;
end
if(r < 1 | r > Nr)
  fprintf('ERROR: r=%d, out of range [1 %d]\n',r,Nr);
  return;
end
if(c < 1 | c > Nc)
  fprintf('ERROR: c=%d, out of range [1 %d]\n',c,Nc);
  return;
end

i = sub2ind([Ns Nr Nc],s,r,c);

V = reshape(V,[Nv Nev]);
Vi = V(i,:);

Fcrs = Vi*V'; %'
Fcrs(i) = 1 - Fcrs(i);

Fcrs = reshape(Fcrs,[Ns Nr Nc]);

return;
