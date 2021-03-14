function [s,fa] = ssbloch(tr,te,fa,t1,t2s,pd)
% [s fa] = ssbloch(tr,te,fa,t1,t2s,<pd>)
%
% Steady-state Bloch Equation


% tr = repetition time
% te = echo time
% fa = flip angle (radians)
% t1 = T1
% t2s = T2 star
% pd = proton density (default = 1)
%
% Time units don't matter as long as they are consitent
%
% if fa=[], set to ernst angle:
%   fa = acos(exp(-TR/T1))
%
% From: Wansapura, et al, J MAG RESE IMG 9:531 538 (1999)
%  At 3T, 
%  Gray:  T1 = 1331ms, T2* = 42ms (occipital) - 52ms (frontal)
%  White: T1 =  832ms, T2* = 48ms (occipital) - 44ms (frontal)
%  CSF:   T1 = 4163ms (Chen Proc SMRI, 2001), T2=503
%  Caud:  T1 = 1271ms                  (Chen Proc SMRI, 2001)
%   
% Zhang Magnetic Resonance in Medicine 70:1082–1086 (2013)
% T1 for blood (Males and Females)
% 1.5T   1531 1429
% 3.0T   1681 1618
% 7.0T   2163 2012
% 
% Yuval Zur, Saul Stokar, Peter Bendel. An analysis of fast imaging
% sequences with steady-state transverse magnetization refocusing.
% MRM, 6:2 175-193, 1988.
%
% Hai-Ling Margaret Cheng, Graham A Wright. Rapid high-resolution T1
% mapping by variable flip angles: Accurate and precise measurements
% in the presence of radiofrequency field inhomogeneity. MRM, Volume
% 55 Issue 3, Pages 566 - 574, 2006.

%
% ssbloch.m
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

% This might also be a useful ref:
% Steen, et al, Precise and accurate measurement of proton T1 in human
% brain in vivo: Validation and preliminary clinical application
% Journal of Magnetic Resonance Imaging Volume 4 Issue 5, Pages 681 -
% 691.

% Rooney, 2007, Magnetic field and tissue dependencies of human brain
% longitudinal1H2O relaxation in vivo Magnetic Resonance in Medicine,
% Volume 57, Issue 2, 308-318. CSF T1=4300

% Voxel-based analysis of R2* maps in the healthy human brain Journal
% of Magnetic Resonance Imaging, Volume 26, Issue 6, 2007, First Page
% 1413 Peran.

% RapidT1 mapping using multislice echo planar imaging Magnetic
% Resonance in Medicine, Volume 45, Issue 4, 2001, First Page 630
% Clare, Stuart; Jezzard, Peter. T1 of CSF = 3700 +/- 500.

% This is something from Gary G. Mildly related.  Using IR, collect
% data with TR >> expected_T1, e.g. 5s, and with TIs of 100 200 500
% 1000 3000 ms.  Fit to
% 
% S(n) = So*(1-2*(1+eps)*exp(-TI(n)/T1) for every voxel,
% 
% where So is proton density, and eps is an error term that accounts
% for inaccuracies in flip angle.  I have a program that gives me the
% three images (T1, So, eps).  eps is usually < 0.05, but depends on
% coil and Bo.


s = [];
if(nargin < 5 | nargin > 6)
  fprintf('[s fa]= ssbloch(tr,te,fa,t1,t2s,<p>)\n');
  return;
end

if(~exist('pd','var')) pd = 1; end
if(isempty(fa)) fa = acos(exp(-tr./t1)); end

etrt1 = exp(-tr./t1);
s = pd .* sin(fa) .* (1-etrt1) .* exp(-te./t2s) ./ (1-cos(fa).*etrt1);

return;






