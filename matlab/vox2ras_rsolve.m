function [M_v1, M_v2, M_Ru1, M_Ru2] = vox2ras_rsolve(Vc_C, inPlaneRotation)
%%
%% NAME
%%
%%     vox2ras_rsolve.m (vox2ras_r{otation matrix}solve) 
%%
%% AUTHOR, VERSION
%%
%%	Rudolph Pienaar, $Id: vox2ras_rsolve.m,v 1.2 2004/05/26 20:16:28 rudolph Exp $
%%
%% SYNOPSIS
%%
%%     [M_v1 M_v2 M_Ru1 M_Ru2] = vox2ras_rsolve(Vc_C, inPlaneRotation)
%%
%% ARGUMENTS
%%
%%      Vc_C		in      column vector defining a direction cosine for
%%					a volume
%%	inPlaneRotation	in	scalar float defining the in plane rotation as
%%					read from the data meas.asc file
%%      M_v1		out     1st candidate vox2ras matrix (no rotation)
%%	M_v2		out	2nd candidate vox2ras matrix (no rotation)
%%	M_Ru1		out	v1 rotated by f(inPlaneRotation)
%%	M_Ru2		out	v2 rotated by f(inPlaneRotation)
%%
%%		where f() is a linear function on inPlaneRotation
%%
%% DESCRIPTION
%%
%%	"vox2ras_rsolve" attempts to find a candidate vox2ras rotation matrix that 
%%	contains the vector C = [ci cj ck]' such that:
%%	
%%				ai	bi	ci 	0
%%		vox2ras	=	-1	bj	cj	0
%%				ak	-1	ck	0
%%				 0	 0	 0	1
%%
%%	Two candidate matrices are returned due to the quadratic nature of the
%%	solutions.
%%
%% PRECONDITIONS
%%
%%	o The vector C is read from a Siemens meas.asc file such that
%%		ci	= sSliceArray.asSlice[0].sNormal.dSag
%%		cj	= sSliceArray.asSlice[0].sNormal.dCor
%%		ck	= sSliceArray.asSlice[0].sNormal.dTra
%%
%% POSTCONDITIONS
%%
%%	o All returned matrices are 4x4.
%%	o Only the rotations of the vox2ras matrix are determined by this function.
%%		The center of k-space is not determined.
%%	o M_v1 M_v2: 2 vox2ras matrices corresponding to two quadratic solutions.
%%	o M_Ru1 M_Ru2: v1 v2 rotated by f(inPlaneRotation).
%%
%% SEE ALSO
%%
%%	vox2ras_ksolve	- determine the k-space col in RAS of a vox2ras matrix
%%	vox2ras_dfmeas	- main function: determines the vox2ras matrix from a
%%			  Siemens meas.asc file.
%% 
%% HISTORY
%%
%% 18 May 2004
%% o Initial design and coding.
%%
%% 26 May 2004
%% o Touching up... 
%% o 4x4 return sizes fixed.
%%

ci	= Vc_C(1);
cj	= Vc_C(2);
ck	= Vc_C(3);
vox2ras	= zeros(3, 3);
v1	= zeros(3, 3);
v2	= zeros(3, 3);

%% We know that each column vector in vox2ras = [A B C] is orthogonal, i.e. the
%% dot products of each column vector is zero (A.B = B.C = A.C = 0). We also
%% know that C = B X A. Thus, we can expand the cross and dot products and solve
%% for ak first, arriving at the following quadratic equation:
%%
%%	(ci^2 + ck^2)ak^2 - (2cj ck + cicjck)ak + (ci+1)(cj^2 + ci^2) = 0

a	= ci^2 + ck^2;
b	= -(2*cj*ck + ci*cj*ck);
c	= (ci+1)*(cj^2 + ci^2);

ak1	= (-b + sqrt(b^2-4*a*c))/2/a;
ak2	= (-b - sqrt(b^2-4*a*c))/2/a;

ai1	= (cj - ak1*ck)/ci;
ai2	= (cj - ak2*ck)/ci;


bj1	= (ci + 1) / ak1;
bi1	= (ck - cj*bj1) / ci;
bj2	= (ci + 1) / ak2;
bi2	= (ck - cj*bj2) / ci;

M3_v1 = [
	ai1	bi1	ci
	-1	bj1	cj
	ak1	-1	ck
];

for col=1:3,
	M_v1(:,col) = M_v1(:,col) ./ norm(M_v1(:,col));
end

M3_v2 = [
	ai2	bi2	ci
	-1	bj2	cj
	ak2	-1	ck
];

for col=1:3,
	M_v2(:,col) = M_v2(:,col) ./ norm(M_v2(:,col));
end

%% The above calculated rotation matrices define an (x,y) plane
%%	given by the first two column vectors. Since we have
%%	"hardcoded" two components of these vectors, we have
%%	fixed them along a default orientation. This default
%%	orientation needs to be rotated by an additional theta_f
%%	radians in order to arrive at the final vox2ras matrix.
%%	By studying several existing vox2ras matrices and their
%%	meas.asc InPlaneRotation values, a simple linear relationship
%%	between this theta_f and the InPlaneRotation value was
%%	found:
%%
%%		theta_f	= m * InPlaneRotation + b
%%
%%	where
%%		m = -1.0025
%%		b = -0.5188 (roughly equal to pi/6 = )
%%

m		= -1.0025;
b		= -0.5188;
theta_f		= m*inPlaneRotation + b;
%theta_f		= inPlaneRotation;

M3_Mu	= [	 cos(theta_f)	 sin(theta_f)	0
		-sin(theta_f)	 cos(theta_f)	0
		 	0		0	1];
	
M3_Ru1	= M3_v1 * M3_Mu;
M3_Ru2	= M3_v2 * M3_Mu;

M_v1	= eye(4);	M_v1(1:3, 1:3)	= M3_v1;
M_v2	= eye(4);	M_v2(1:3, 1:3)	= M3_v2;
M_Ru1	= eye(4);	M_Ru1(1:3, 1:3)	= M3_Ru1;
M_Ru2	= eye(4);	M_Ru2(1:3, 1:3)	= M3_Ru2;

%% All done!
