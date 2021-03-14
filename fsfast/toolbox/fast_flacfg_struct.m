function flacfg = fast_flacfg_struct
% Create an FS-FAST first-level analysis structure with all 
% the fields. 


%
% fast_flacfg_struct.m
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

flacfg.version = 2;
flacfg.flaname = '';
flacfg.TR = [];
flacfg.fxlist = [];          % Effects Models, see fast_fxcfg_struct.m
flacfg.inorm = 0;       % Intensity Normalization Target
flacfg.nskip = [];
flacfg.slicetiming = '';
flacfg.noisemodel = [];

flacfg.evschfname = '';  % Event Schedule File Name, Relative
flacfg.runlistfile = ''; % Run List File, Relative
flacfg.runweightfile = ''; % Run Weight File, Relative
flacfg.fsd = '';
flacfg.funcstem = '';   
flacfg.maskstem = '';  % Relative to FSD/masks

flacfg.useevschweight = []; % Use Event Weights, if present in EvSch
flacfg.usetpexclude = [];   % Use TP Exclude, if present in EvSch

% These variable do not really describe the analysis, but are 
% a convenient way for for keeping state
flacfg.sesspath = ''; % Absolute session path
flacfg.sesscfg = [];  % See fast_sesscfg_struct
flacfg.nthfx = [];
flacfg.nthrun = [];
flacfg.tDelay = [];

return;
