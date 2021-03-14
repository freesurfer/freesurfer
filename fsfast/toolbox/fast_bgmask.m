function [bgm, fstdbg, fstdfg, fstdfgexp, bgthresh, fgmn, bgmB] = fast_bgmask(fmn,fstd,noepi)
% [bgm, fstdbg, fstdfg, fstdfgexp, bgthresh, fgmn, bgmB] = fast_bgmask(fmn,fstd,<noepi>)
%
% White noise floor in the foreground should then be
%     fstdfgexp = fstdbg/sqrt(2-pi/2);
%
% bgmB is the background mask as determined by Stoecker's method.
%
%
%


%
% fast_bgmask.m
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

bgm = [];
fstdbg = [];
fstdfg = [];
fstdfgexp = [];
bgthresh = [];
fgmn = []; 
bgmB = [];
if(nargin < 2 | nargin > 3)
  fprintf('[bgm fstdbg fstdfg fstdfgexp bgthresh] = fast_bgmask(fmn,fstd,<noepi>)\n');
  return;
end

if(nargin == 2) noepi = 0; end

indz = find(fstd==0);
gmean = mean(fmn(:));

% First, segment the foreground
fgm = fmn > gmean; % Threshold mean image
fgm = fast_dilate(fgm,3,1); % Erode by 3
fgm(:,:,1) = 0;   % Exclude first slice
fgm(:,:,end) = 0; % Exclude last slice
indfg = find(fgm);
fstdfg = sqrt(mean(fstd(indfg).^2));
fgmn = mean(fmn(indfg));

% Use Foreground std to set the background thresh
bgthresh = 4*fstdfg*sqrt(pi/2);

% Now create a new foreground mask
fgm2 = fmn > bgthresh;
fgm2 = fast_dilate(fgm2,5); % Dilate by 5

if(~noepi)
  % Ghost of the fg mask
  fggm = fast_ghostmask(fgm2); 
  % Create a mask of everything not wanted in background
  m = fgm2 | fggm; % Foreground and ghost
else
  m = fgm2;
end


m(1,:,:)   = 1;  % Set first row to 1
m(end,:,:) = 1;  % Set last  row to 1
m(:,1,:)   = 1;  % Set first col to 1
m(:,end,:) = 1;  % Set last  col to 1
m(indz)    = 1;  % Excl vox with 0 std

% Now the background mask is everything else
bgm = ~m; 

indbg = find(bgm);
nbgm = length(indbg);
if(nbgm == 0)
  fprintf('ERROR: could not segment background\n');
  keyboard
  return;
end

fstdbg = sqrt(mean(fstd(indbg).^2));
fstdfgexp = fstdbg/sqrt(2-pi/2);

% Just for fun, use the Stoecker method
bgmB = zeros(size(m));
dr = floor(size(bgmB,1)/8);
dc = floor(size(bgmB,2)/8);
ds = floor(size(bgmB,3)/8);
bgmB(1:dr,         1:dc, 1:ds) = 1;
bgmB(end-dr+1:end, 1:dc, 1:ds) = 1;
bgmB(1:dr,         end-dc+1:end, 1:ds) = 1;
bgmB(end-dr+1:end, end-dr+1:end, 1:ds) = 1;
bgmB(1:dr,         1:dc, end-ds+1:end) = 1;
bgmB(end-dr+1:end, 1:dc, end-ds+1:end) = 1;
bgmB(1:dr,         end-dc+1:end, end-ds+1:end) = 1;
bgmB(end-dr+1:end, end-dr+1:end, end-ds+1:end) = 1;
% Remove all edges and vox with 0 std
bgmB(1,:,:)   = 0;
bgmB(end,:,:) = 0;
bgmB(:,1,:)   = 0;
bgmB(:,end,:) = 0;
bgmB(indz)    = 0;

return;

