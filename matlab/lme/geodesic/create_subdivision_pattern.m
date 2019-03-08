%regular subdivision pattern of a triangle
%used in drawing approximate equidistant lines in example5.m
%Copyright (c) 2007 Danil Kirsanov
function [weights,tri] = create_subdivision_pattern(level)    %"level" is the number of additional vertices per edge

step = 1/(level + 1);
N = level + 2;                      %total number of points per edge

N_p = (N^2 + N)/2;                  %number of vertices in the regular grid
N_t = (N-1)^2;                      %two triangles per square
weights = zeros(N_p,3);
tri = zeros(N_t,3);

n_p = 1;
n_t = 1;
for i=1:N;
    K = N-i+1;
    for j=1:K
        c = step*[j-1, i-1];
        weights(n_p,:) = [c, 1-sum(c)];       
       
        if i<N & j<K
            tri(n_t,:) = [n_p, n_p + 1, n_p + K];
            n_t = n_t + 1;
            if j<K-1
                tri(n_t,:) = [n_p + 1, n_p + K, n_p + K + 1];
                n_t = n_t + 1;
            end;
        end;
        
        n_p = n_p+1;        
    end;
end;