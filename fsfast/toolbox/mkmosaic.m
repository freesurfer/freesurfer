% mkmosiac.m
%
% Makes a NxMxnTP mosaic from input slice files. If the
% slice files have nTP layers, then the mosaic will have nTP layers.
% Used by the csh script mkmosaic.
%
% $Id: mkmosaic.m,v 1.1 2003/03/04 20:47:41 greve Exp $


fprintf('\n');

if( ~exist('SliceFile') )
  fprintf('ERROR: variable SliceFile does not exist\n');
  qoe;
  return;
end
if( ~exist('MosaicFile') )
  fprintf('ERROR: variable MosaicFile does not exist\n');
  qoe;
  return;
end
if( ~exist('HistEQ') )
  fprintf('ERROR: variable HistEQ does not exist\n');
  qoe;
  return;
end
if( ~exist('HistEQThr') )
  fprintf('ERROR: variable HistEQThr does not exist\n');
  qoe;
  return;
end
if( ~exist('MosaicCols') )
  fprintf('ERROR: variable MosaicCols does not exist\n');
  qoe;
  return;
end

clear y;

nSlices = size(SliceFile,1);

for n = 1:nSlices,

  s = deblank(SliceFile(n,:));
  fprintf('Loading %s ...\n',s);
  slice = fmri_ldbfile(s);
  if(isempty(slice))
    fprintf('Error loading %s\n',s);
    qoe;
    return;
  end

  y(:,:,:,n) = slice;

end

fprintf('Making mosiac with %d columns\n',MosaicCols);
ymos = mstk2mos(y,MosaicCols);

if(HistEQ)
  for n = 1:size(ymos,3),
    fprintf('Equalizing Mosaic Slice %2d with threshold %g ...\n',n,HistEQThr);
    [ymos(:,:,n) xthresh nclip psqueeze] = drsqueeze(ymos(:,:,n),HistEQThr);
    fprintf('  xthresh  = %g\n',xthresh);
    fprintf('  nclip    = %d / %d\n',nclip,prod(size(ymos(:,:,n))));
    fprintf('  psqueeze = %g\n',psqueeze);
  end
end

fprintf('Saving mosiac in %s\n',MosaicFile);
fmri_svbfile(ymos,deblank(MosaicFile));

