function samseg_writeOutFreeSurferSeg(imageFileName,transformedTemplateFileName,croppedBuffer,savePath)
% This is Oula's code renamed. It basically uncrops to be the resolution of the input.
%
% $Id: samseg_writeOutFreeSurferSeg.m,v 1.1 2017/01/26 00:21:49 greve Exp $
%

%First get the cropping indices out
[origInd, croppedInd, croppedDIM, origDIM] = kvlGetCroppedRegion(imageFileName,transformedTemplateFileName);

%Next create a buffer of the original image size
origBuf = zeros(origDIM,'single');

%Then set the cropped buffer to the right place
%croppedInd = [size(croppedBuffer,1)-croppedDIM(1) size(croppedBuffer,2)-croppedDIM(2) size(croppedBuffer,3)-croppedDIM(3)];

origBuf(origInd(1)+1:origInd(1)+croppedDIM(1),origInd(2)+1:origInd(2)+croppedDIM(2),origInd(3)+1:origInd(3)+croppedDIM(3))=...
croppedBuffer(croppedInd(1)+1:croppedInd(1)+croppedDIM(1),croppedInd(2)+1:croppedInd(2)+croppedDIM(2),croppedInd(3)+1:croppedInd(3)+croppedDIM(3));

%Get the header information
im = kvlReadImage(imageFileName);
kvlSetImageBuffer(im,origBuf);

kvlWriteImage(im,[savePath '/segSubSpace.mgz']);

return

