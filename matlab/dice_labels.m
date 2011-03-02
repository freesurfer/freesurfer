function dice = dice_labels(lname1, lname2)
% d = dice_labels(lname1, lname2)
%
% computes Dice coefficient of the two given label files.
% 

%
% dice_labels.m
%
% Original Author: Nick Schmansky
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
%    $Revision: 1.2 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

label1 = read_label('',lname1);
label2 = read_label('',lname2);

label1Size = size(label1,1);
label2Size = size(label2,1);

hits = 0;
for i=1:label1Size
   x1 = label1(i,2);
   y1 = label1(i,3);
   z1 = label1(i,4);
   for j=1:label2Size
     x2 = label2(j,2);
     y2 = label2(j,3);
     z2 = label2(j,4);
     if ((x1==x2) && (y1==y2) && (z1==z2))
       hits = hits+1;
       break;
     end
   end
end

dice = (2 * hits) / (label1Size + label2Size);


