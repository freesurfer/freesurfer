function [ cost reshapedGradient ] = kvlEvaluateMeshPositionInVectorFormat( reshapedPosition, mesh, image, transform, means, variances )
%
% Simple function that evaluates the cost and gradient of the location of the nodes of 
% tetrahedral mesh-based atlas.
% 
% The only thing this really does is convert position and gradients from [ numberOfNodes 3 ] format
% into vectorized format, so that we can use general-purpose multidimensional numerical optimizers. 
%
kvlSetMeshNodePositions( mesh, reshape( reshapedPosition, [ length( reshapedPosition )/3 3 ] ) );
[ cost gradient ] = kvlEvaluateMeshPosition( mesh, image, transform, means, variances );
reshapedGradient = gradient(:);

