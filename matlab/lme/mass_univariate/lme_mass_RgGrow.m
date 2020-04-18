function [Regions,RgMeans] = lme_mass_RgGrow(SphSurf,Re,Theta,maskvtx,nst,prc)
% [Regions,RgMeans] = lme_mass_RgGrow(SphSurf,Re,Theta,maskvtx,nst,prc)
% 
% This function implements a region growing algorithm along the spherical
% surface to find homogeneous regions comprising locations with similar  
% covariance components.
%
% Input
% SphSurf: Spherical surface, which represent a spherical coordinate system
% for the surface where the analysis is to be done. It is a structure in
% fs_read_surf format. SphSurf.tri = t x 3 matrix of triangle indices,
% 1-based, t=#triangles and SphSurf.coord = 3 x nv matrix of coordinates,
% nv=#vertices.
% Re: Matrix of residual errors at each location (nmxnv, nm total # of
% maps, nv #vertices.
% Theta: Matrix whose colums are estimators of the covariance components at.
% each location.
% maskvtx: Mask's vertices (1-based). Default [] (all vertices included).
% nst: Number of standard deviations of Theta inside a region for it to
% be considered as homogeneous. Default 2.
% prc: Percent of vertices inside an homogeneous region that are allowed to
% not meet the homogenity criteria. This parameter makes the algorithm
% robust to noise and outliers. Default 95%.
%
% Output
% Regions: 1 x nv segmentation vector containing a region number assigned 
% to each vertice along the surface.
% RgMeans: Matrix assigning the mean of Theta within a region to each vertex
% belonging to that region.
%
% Original Author: Jorge Luis Bernal Rusiel 
% References: 
% Gonzalez, R.C., Woods, R.E. (2002). Digital Image Processing. 2nd Edition, 
% New Jersey: Prentice Hall, Inc.
%
tic;
if nargin < 3
    error('Too few inputs');
elseif nargin < 6
        prc = 95;
        if nargin < 5
            nst = 2;
            if nargin < 4
                maskvtx = [];
            end
        end;
end;
nv = size(Theta,2);
if isempty(maskvtx)
    maskvtx = 1:nv;
end;
AdjM = AdjMtx(SphSurf,maskvtx);
display(' ');
display('Starting region growing algorithm');
[phi,theta] = sphCoord(SphSurf.coord');
sphcoord = [phi,theta]';
Regions = vtxwGrowing(Re,Theta,sphcoord,AdjM,maskvtx,nst,prc);
[RgMeans,nRg] = lme_mass_RgMean(Regions,Theta);
display(' ');
display([num2str(nRg-1) ' homogeneous regions were obtained.']);
et = toc;
display(['Elapsed time is ' num2str(et/60) ' minutes.']);
end






function Rgs = vtxwGrowing(Re,Params,coord,Adj,maskvtx,nst,prc)
display('Computing seeds ...');
[Rgs1,Rgseed,Rgstd]= seeds(Re(:,maskvtx),Params(:,maskvtx),coord(:,maskvtx),nst,prc);
Rgseed = maskvtx(Rgseed);
ns = length(Rgseed);
display([num2str(ns) ' seeds were computed after splitting the surface']);
[np,nv] = size(Params);
Rgs = zeros(1,nv,'uint32');
Rgs(maskvtx) = Rgs1;
notIncl = true(1,nv);
notIncl(Rgs == 0) = false;
notIncl(Rgseed) = false;
maxRgsz = 90; %Do not form regions with size greater than maxRgsz 
Rgvtxs = zeros(ns,maxRgsz,'uint32');
Rgvtxs(:,1) = Rgseed;
Rgvtxs_prvtxind = zeros(1,ns,'uint32');
Rgvtxs_lastind = ones(1,ns,'uint32');
Rgind = 1:ns;
thr = nst * Rgstd;
growing_ns = ns;
while (growing_ns > 0)
    i = 1;
    display(' ');
    display([num2str(growing_ns) ' seeds for growing']);
    display(['Current maximum region size ' num2str(max(Rgvtxs_lastind)) ' vertices']);
    while (i <= growing_ns)
        Rgvtxs_prvtxind(Rgind(i)) = Rgvtxs_prvtxind(Rgind(i)) + 1;
        prvtx = Rgvtxs(Rgind(i),Rgvtxs_prvtxind(Rgind(i)));
        RgParams = Params(:,Rgvtxs(Rgind(i),1:Rgvtxs_lastind(Rgind(i))));
        mRgParams = mean(RgParams,2);
        Rgsz = Rgvtxs_lastind(Rgind(i));
        %processing the current vertex prvtx
        adjvtxs = Adj(prvtx,Adj(prvtx,:) > 0);
        %try to include only the adjacent vertices not included in other
        %regions or in this region itself
        adjvtxs = adjvtxs(notIncl(adjvtxs));
        %try to include first the adjacent vertices with the param closest
        %to the region mean param.
        [~,loc] = min(sum(abs(Params(:,adjvtxs) - kron(ones(1,length(adjvtxs)),mRgParams)),1));
        NewRgParams = [RgParams Params(:,adjvtxs(loc))];
        mNewRgParams = mean(NewRgParams,2);
        nRgnv = size(NewRgParams,2);
        DistNewRgParams = abs(NewRgParams-kron(ones(1,nRgnv),mNewRgParams));
        %an adjacent vertex is included in the region if the difference
        %between its params and the region mean params is bellow the
        %threshold thr or if the region still contains a percent
        %prc of vertices bellow the threshold.
        nadjv = length(adjvtxs); j = 1;
        while (Rgsz < maxRgsz) && (j <= nadjv) && ...
                sum(sum(DistNewRgParams <= kron(ones(1,nRgnv),thr(:,Rgind(i))),1) == np) >= prc*(nRgnv)/100
            %In addition the vertex's residuals must be correlated above 0.5
            %with any other vertex's residual in the region.
            if min(min(corrcoef([Re(:,Rgvtxs(Rgind(i),1:Rgvtxs_lastind(Rgind(i)))),Re(:,adjvtxs(loc))]))) >= 0.5
                Rgvtxs_lastind(Rgind(i)) = Rgvtxs_lastind(Rgind(i)) + 1;
                Rgsz = Rgsz + 1;
                Rgvtxs(Rgind(i),Rgvtxs_lastind(Rgind(i))) = adjvtxs(loc);
                notIncl(adjvtxs(loc)) = false;
                RgParams = Params(:,Rgvtxs(Rgind(i),1:Rgvtxs_lastind(Rgind(i))));
                mRgParams = mean(RgParams,2);
            end;
            adjvtxs = [adjvtxs(1:loc-1),adjvtxs(loc+1:end)];
            [~,loc] = min(sum(abs(Params(:,adjvtxs) - kron(ones(1,length(adjvtxs)),mRgParams)),1));
            NewRgParams = [RgParams Params(:,adjvtxs(loc))];
            mNewRgParams = mean(NewRgParams,2);
            nRgnv = size(NewRgParams,2);
            DistNewRgParams = abs(NewRgParams-kron(ones(1,nRgnv),mNewRgParams));
            j = j + 1;
        end;
        if (Rgvtxs_prvtxind(Rgind(i)) == Rgvtxs_lastind(Rgind(i))) || (Rgvtxs_lastind(Rgind(i)) == maxRgsz)
            Rgind = [Rgind(1:i-1) Rgind(i+1:end)];
            growing_ns = growing_ns - 1;
        else
            i = i + 1;
        end;
    end;
end;
for i=1:ns
    Rgs(Rgvtxs(i,1:Rgvtxs_lastind(i))) = i;
end;
nni = sum(notIncl);
notInclv = find(notIncl == 1);
Rgs(notInclv) = ns+1:ns+nni;
%try to add holes to their most similar adjacent homogeneous region
display(' ');
display([num2str(nni) ' holes (unassigned vertices). Trying to form new'...
        ' regions among them or add them to their most similar adjacent'...
        ' region ...']);
for i=1:nni
    hadjvtxs = Adj(notInclv(i),Adj(notInclv(i),:) > 0);
    if ~isempty(hadjvtxs)
        closestdistRParams = sum(abs(Params(:,notInclv(i)) - mean(Params(:,Rgs == Rgs(hadjvtxs(1))),2)));
        nadjvtxs = length(hadjvtxs);
        m = 1;
        j = 2;
        while j <= nadjvtxs
            distRParams = sum(abs(Params(:,notInclv(i)) - mean(Params(:,Rgs == Rgs(hadjvtxs(j))),2)));
            if distRParams < closestdistRParams
                closestdistRParams = distRParams;
                m = j;
            end;
            j = j+1;
        end;
        [Rgs,mrg] = merge(Rgs,Re,Params,Rgs(notInclv(i)),Rgs(hadjvtxs(m)),nst,prc);
        if mrg
            notIncl(notInclv(i)) = 0;
        end;
    end;
end;
nni = sum(notIncl);
if nni > 0
    display([num2str(nni) ' holes (unassigned vertices) remained. Will be'...
        ' treated as independent regions of size 1.']);
else
    display(' No holes (unassigned vertices) remained. ');
end;

end






function [Rgs,Rgseed,Rgstd] = seeds(Re,Params,coord,nst,prc)
[np,nv] = size(Params);
Rgs = ones(1,nv,'uint32');
spl = true;
while spl 
    splRgs = unique(Rgs);
    ntospl = length(splRgs);
    spltf = zeros(1,ntospl);
    nRgs = ntospl;
    for i=1:ntospl
        [Rgs,spltf(i),nnRgs] = split(Rgs,Re,Params,coord,splRgs(i),nst,prc);
        nRgs = nRgs + nnRgs;
    end;
    if sum(spltf) == 0
        spl = false;
    end;
end;
    splRgs = unique(Rgs);
    nRgs = length(splRgs);
    Rgseed = zeros(1,nRgs);
    Rgstd = zeros(np,nRgs);
    for i=1:nRgs
        Rgv = find(Rgs == splRgs(i));
        RgParams = Params(:,Rgv);
        Rgstd(:,i) = std(RgParams,0,2);
        nv = length(Rgv);
        spp = sum(coord(:,Rgv),2)/nv;
        [~,loc] = min(sum(abs(coord(:,Rgv) - kron(ones(1,nv),spp)),1));
        Rgseed(i) = Rgv(loc);
    end;
end






function [Rgs,tf,nRnv] = split(Rgs,Re,Params,coord,Rnum,nst,prc)
if homogenity(Re,Params,Rgs,Rnum,nst,prc)
    nRnv = 0;
    tf = false;
else
    %split the non-homogeneous region
    Rvts = find(Rgs == Rnum);
    nv = length(Rvts);
    spp = sum(coord(:,Rvts),2)/nv;
    nR = max(Rgs);
    for i=1:nv
        if (coord(1,Rvts(i)) <= spp(1)) && (coord(2,Rvts(i)) <= spp(2))
            Rgs(Rvts(i)) = 1;
        elseif (coord(1,Rvts(i)) <= spp(1)) && (coord(2,Rvts(i)) > spp(2))
            Rgs(Rvts(i)) = 2;
        elseif (coord(1,Rvts(i)) > spp(1)) && (coord(2,Rvts(i)) <= spp(2))
            Rgs(Rvts(i)) = 3;
        elseif (coord(1,Rvts(i)) > spp(1)) && (coord(2,Rvts(i)) > spp(2))
            Rgs(Rvts(i)) = 4;
        end;
    end
    Rnv = unique(Rgs(Rvts));
    nRnv = length(Rnv) - 1;
    for i=1:nRnv
        Rgs(Rvts(Rgs(Rvts) == Rnv(i+1))) = nR+i;
    end;
    Rgs(Rvts(Rgs(Rvts) == Rnv(1))) = Rnum;
    mrgRgs = [Rnum nR+1:nR+nRnv];
    %merge homogeneous subregions
    nmtomrg = length(mrgRgs);
    i = 1; 
    while i <= nmtomrg
        j = i + 1; mrg = false;
        while (~mrg) && (j <= nmtomrg)
            [Rgs,mrg] = merge(Rgs,Re,Params,mrgRgs(i),mrgRgs(j),nst,prc);
            j = j + 1;
        end;
        if mrg 
            mrgRgs = [mrgRgs(1:j-2) mrgRgs(j:end)];
            nmtomrg = nmtomrg - 1;
        else
            i = i + 1;
        end;
    end;
    tf = true;
end;
end





function tf = homogenity(Re,Params,Rgs,Rn,nst,prc)
tf = false;
maxRgsz = 90; 
%Regions with size greater than maxRgsz are not considered homogeneous
Rgv = find(Rgs == Rn);
if length(Rgv) <= maxRgsz
    RgParams = Params(:,Rgv);
    mRgParams = mean(RgParams,2);
    stRgParams = std(RgParams,0,2);
    thr = nst*stRgParams;
    [np,nv] = size(RgParams);
    DistRgParams = abs(RgParams-kron(ones(1,nv),mRgParams));
    if (sum(sum(DistRgParams <= kron(ones(1,nv),thr),1) == np) >= prc*nv/100) &&...
            min(min(corrcoef(Re(:,Rgv)))) >= 0.5
        tf = true;
    end;
end;
end





function [Rgs,tf] = merge(Rgs,Re,Params,Rn1,Rn2,nst,prc)
if Rn1 <= Rn2
    Rkn = Rn1;
    Ren = Rn2;
else
    Rkn = Rn2;
    Ren = Rn1;
end;
auxRgs = Rgs;
auxRgs(auxRgs == Ren) = Rkn;
if homogenity(Re,Params,auxRgs,Rkn,nst,prc)
    Rgs = auxRgs;
    tf = true;
else
    tf = false;
end;
end
