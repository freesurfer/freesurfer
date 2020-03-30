function [sX,sY,ni,ssID] = sortData(X,tcol,Y,sID)
% [sX,sy,ni,ssID] = sortData(X,tcol,y,sID)
%
% This function sorts the design matrix X and the data Y according to 
% subject's ID (sID) and according to time within each subject. 
%
% Input
% X: Design Matrix (nmxp, nm total # of measurements, p #of fixed effects). 
% tcol: Column number in X for the time covariate. 
% Y: Data matrix (nmxnv, nv #locations) whos colums can represent data 
% vectors for each voxel/vertex.
% sID: Vector or cell array of subjects'ID (nmx1). This is a vector 
% indicating a subject ID for each row of X and Y. 
%
% Output
% sX: Ordered design matrix.
% sY: Ordered data matrix.
% ni: Vector (mx1, m # of subjects in the study) whose entries are the 
% number of repeated measures for each subject (ordered according to sX and
% sY).
% ssID: Vector of ordered subjects'ID. 
%
% Original Author: Jorge Luis Bernal Rusiel 
%  
if nargin < 4 
    error('Too few inputs');   
end;
[nm,p] = size(X);
if (tcol<1) || (tcol>p)
      error(['The time column input (tcol) must be greater than 0 and less'...
         ' than the number of colums in X (' num2str(p) ')']);   
end; 
nY = size(Y,1);
if (nY~=nm) || (length(sID)~=nm)
      error(['The number of rows in the inputs Y and sID must be the same'...
         ' as the number of rows in X']);   
end;  
[ssID,sID_ind] = sort(sID);
count = 1;
m = length(unique(ssID));
ni = ones(m,1);
for i=2:nm
    if strcmp(char(ssID(i-1)),char(ssID(i)))
        ni(count) = ni(count) + 1;
    else
      count = count + 1;  
    end;     
end;
time = X(sID_ind,tcol);
sX = X;
sY = Y;
count = 0;
for i=1:m
   [~,time_ind] = sort(time(count+1:count+ni(i)));
   sX(count+1:count+ni(i),:) = X(sID_ind(count+time_ind),:);
   sY(count+1:count+ni(i),:) = Y(sID_ind(count+time_ind),:);
   count = count + ni(i);
end;
  
