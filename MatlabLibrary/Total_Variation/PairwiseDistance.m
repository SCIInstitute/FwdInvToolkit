function [D] = PairwiseDistance(varargin)
% PairwiseDistance This function calculates pairwise distances (2norm) between points
% in X and Y or, if this function is called with one argument, between
% points in X only.
%   INPUT ARGS
%       varargin{1}: matrix with points as columns
%       varargin{2}: matrix with points as columns
%   OUTPUT ARGS
%       D: distance matrix. D(i,j) corresponts to the distance between the
%       ith point in varargin{1} and the jth point in varargin{2}. I f
%       there is only one argument, D(i,j) is the distance between ith
%       point in varargin{1} and jth point in varargin{1}.
    
    if(nargin ==1)
        %% Calculates the pairwise distances between points in a matrix X
        % dist(xi,xj) = (xi - xj)^2 = ||xi||^2 + ||xj||^2 - 2*xi'*xj
        % //2norm
        
        K = varargin{1}'*varargin{1};
        vD = diag(K);
        N = numel(vD);
        A = repmat(vD(:),1,N);
        D = -2*K + A + A';
    
    else
        %% Calculates the pairwise distances between points in matrix X and Y
        % dist(xi,yj) = (xi - yj)^2 = ||xi||^2 + ||yj||^2 - 2*xi'*yj
        % //2norm
        Dy = repmat(sum(varargin{2}.^2,1),size(varargin{1},2),1);
        Dx = repmat(sum(varargin{1}.^2,1)',1,size(varargin{2},2));
        Dxy = varargin{1}'*varargin{2};
        D = Dy + Dx - 2*Dxy;
        
    end
		
end
