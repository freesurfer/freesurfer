function [P,G] = fs_signrank(x_fname,y_fname,varargin)
%
%  Author: Martin Reuter
%
%  Extends Wilcoxon signed rank test for zero median to FS surfaces
%  If y_fname is passed, do a paired test, see fs_ranksum for independent
%  samples.
%
%  Can work in parallel mode (run 'matlabpool OPEN' first)
%
%  x_fname    file name of mgh file (stacked surfaces)
%  y_fname    optional: file name for second file to do a paired test
%
%  fs_signrank(...,'mask',MASK) selects a mask (e.g. cortex.label)
%  
%  fs_signrank(...,'outsig',OUTSIG) specifies file to output p-values
%      output will be in FreeSurfer mgh format, signed - log10(p)
% 
%  fs_signrank(...,'outgamma',OUTGAMMA specifies file to output gamma
%      output will be in FreeSurfer mgh format (see G below)
% 
%  See signrank for more arguments ('method' and 'tail')
%
%  Return (only inside mask):
%
%  G          gamma: median(X)            for single sample
%                    median(X)-median(Y)  for paired samples
%
%  P          p-values of signed rank test
%

% parse and remove mask and outfile params:
oknames = {'mask' 'outsig' 'outgamma' 'tail' 'alpha' 'method'};
dflts   = {'' '' '' 'both' 0.05   '' };
[mask_fname , outsig, outgamma, tail] = internal.stats.parseArgs(oknames,dflts,varargin{:});
% remove mask and outfile:
del = zeros(size(varargin));
for i = 1:length(varargin)
  if (strcmp(varargin{i},'mask'))
    del(i) = 1;
    del(i+1) = 1;
  end
  if strcmp(varargin{i},'outsig')
    del(i) = 1;
    del(i+1) = 1;
  end
  if strcmp(varargin{i},'outgamma')
    del(i) = 1;
    del(i+1) = 1;
  end
end
varargin = varargin(find(del == 0));
%tail = getParamVal(tail,{'both' 'right' 'left'},'''tail''');
if (strcmp(tail,'both'))
  tail = 0;
elseif (strcmp(tail,'left'))
  tail = -1;
elseif (strcmp(tail,'right'))
  tail = 1;
else
  error('tail can only be "both","left" or "right"');
end

% read X from file:
if ( ischar(x_fname) && exist(x_fname,'file') )
  [X,Mx] = fs_read_Y(x_fname);
  nsubj = size(X,1);
  nvert = size(X,2);
else
  error(['Can not load ' x_fname ' as an mgh or mgz file']);
end

% if passed read Y from file:
Y = [];
if (  ~ isempty(y_fname) )
  if ischar(y_fname) && exist(y_fname,'file')
    [Y,My] = fs_read_Y(y_fname);
    if (size(X) ~= size(Y))
      error('Dimensions of X and Y do not agree');
    end
  else    
      error(['Can not load ' y_fname ' as an mgh or mgz file']);
  end
end

% read mask from file:
mask = [1:nvert];
outside = [];
if ( ~ isempty(mask_fname) )
  mask = fs_read_label(mask_fname);
  zz = zeros(1,nvert);
  zz(mask) = ones(1,length(mask));
  outside = find(zz==0);
end

% process
n = length(mask);
Pt = ones(1,n);
if (isempty(Y))
  parfor i = 1:n
    Pt(i) = signrank(X(:,mask(i)),'',varargin{:});
  end
  G = median(X);
else
  parfor i = 1:n
    idx = mask(i);
    Pt(i) = signrank(X(:,idx),Y(:,idx),varargin{:});
  end
  G = median(X)-median(Y);
end
P = ones(1,nvert);
P (mask) = Pt;
G(outside) = zeros(1,length(outside));

% write (on full surface)
if ( ~ isempty(outsig) )
  Pl10 = -log10(P) .* sign(G);
  Mx.volsz(4) = 1;
  fs_write_Y(Pl10,Mx,outsig);
end

if ( ~ isempty(outgamma) )
  Mx.volsz(4) = 1;
  fs_write_Y(G,Mx,outgamma);
end

% return values (only inside mask):
P = Pt;
G = G(mask);
meanG = mean(G);

  
