function [bgm, fstdbg, fstdfg, fstdfgexp, bgthresh] = fast_bgmask(fmn,fstd)
% [bgm, fstdbg, fstdfg, fstdfgexp, bgthresh] = fast_bgmask(fmn,fstd)
%
% White noise floor in the foreground should then be
%     fstdfgexp = fstdbg/sqrt(2-pi/2);
%
% $Id: fast_bgmask.m,v 1.1 2005/03/19 00:22:49 greve Exp $
%

bgm = [];
fstdbg = [];
fstdfg = [];
if(nargin ~= 2)
  fprintf('[bgm fstdbg fstdfg fstdfgexp bgthresh] = fast_bgmask(fmn,fstd)\n');
  return;
end

indz = find(fstd==0);
gmean = mean(fmn(:));

% First, segment the foreground
fgm = fmn > gmean; % Threshold mean image
fgm = fast_dilate(fgm,3,1); % Erode by 3
fgm(:,:,1) = 0;   % Exclude first slice
fgm(:,:,end) = 0; % Exclude last slice
indfg = find(fgm);
fstdfg = sqrt(mean(fstd(indfg).^2));

% Use Foreground std to set the background thresh
bgthresh = 4*fstdfg*sqrt(pi/2);

% Now create a new foreground mask
fgm2 = fmn > bgthresh;
fgm2 = fast_dilate(fgm2,5); % Dilate by 5
fggm = fast_ghostmask(fgm2); % Ghost of the fg mask

% Create a mask of everything not wanted in background
m = fgm2 | fggm; % Foreground and ghost
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
  return;
end

fstdbg = sqrt(mean(fstd(indbg).^2));
fstdfgexp = fstdbg/sqrt(2-pi/2);

return;

