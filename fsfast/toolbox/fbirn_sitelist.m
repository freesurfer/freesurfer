function [sitelist, ind15T, ind3T, ind4T] = fbirn_sitelist
% [sitelist, ind15T, ind3T, ind4T] = fbirn_sitelist
% changed order on 12/19/06


%
% fbirn_sitelist.m
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

sitelist = '';
sitelist = strvcat(sitelist,'bwh');     % 1
sitelist = strvcat(sitelist,'dunc');    % 2
sitelist = strvcat(sitelist,'dunc4t');  % 3
sitelist = strvcat(sitelist,'iowa');    % 4
sitelist = strvcat(sitelist,'mgh');     % 5
sitelist = strvcat(sitelist,'min');     % 6
sitelist = strvcat(sitelist,'nm');      % 7
sitelist = strvcat(sitelist,'stan');    % 8
sitelist = strvcat(sitelist,'uci');     % 9
sitelist = strvcat(sitelist,'ucsd');    % 10
sitelist = strvcat(sitelist,'mgh-te20'); % 11
sitelist = strvcat(sitelist,'mgh-te50'); % 12

ind15T = [2 4 7 9 10];
ind3T = [1 5 6 8 11 12];
ind4T = [3];


return

% This was the original list
% sitelist = '';
% sitelist = strvcat(sitelist,'dunc');    % 1
% sitelist = strvcat(sitelist,'iowa');    % 2
% sitelist = strvcat(sitelist,'nm');      % 3
% sitelist = strvcat(sitelist,'uci');     % 4
% sitelist = strvcat(sitelist,'ucsd');    % 5
% sitelist = strvcat(sitelist,'bwh');     % 6
% sitelist = strvcat(sitelist,'mgh');     % 7
% sitelist = strvcat(sitelist,'min');     % 8
% sitelist = strvcat(sitelist,'stan');    % 9
% sitelist = strvcat(sitelist,'dunc4t');  % 10
