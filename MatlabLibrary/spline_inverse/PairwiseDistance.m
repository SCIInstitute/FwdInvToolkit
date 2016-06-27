function D = PairwiseDistance(varargin)
% PairwiseDistance This function calculates pairwise distances (2norm SQUARED) between points
% in X and Y or, if this function is called with one argument, between
% points in X only.
% Author: Ramon Martinez Orellana
% Updated for large inputs: Burak Erem
% 
% D = PairwiseDistance(A,[B])
%   INPUT ARGS
%       A: matrix with points as columns
%       B: (optional) matrix with points as columns
%   OUTPUT ARGS (NOTE: "distance" here means 2-norm SQUARED of the difference between points)
%       D: distance matrix. D(i,j) corresponts to the distance between the
%       ith point in A and the jth point in B. If
%       there is only one argument, D(i,j) is the distance between ith
%       point in A and jth point in B.
% 
    
    if(nargin == 1)
        % Calculates the pairwise distances squared between columns of a matrix X
        % dist(xi,xj) = (xi - xj)^2 = ||xi||^2 + ||xj||^2 - 2*xi'*xj
        % //2norm squared

        D = distSingle(varargin{1});
    end
    if(nargin == 2)
        % Calculates the pairwise distances between points in matrix X and Y
        % dist(xi,yj) = (xi - yj)^2 = ||xi||^2 + ||yj||^2 - 2*xi'*yj
        % //2norm squared

        D = distDouble(varargin{1},varargin{2});
        
    end
    if(nargin == 3)
        % Block-column-wise computations
        maxcolsize=varargin{3}(1);
        
        N1=size(varargin{1},2);
        N2=size(varargin{2},2);
        
        blocks1=0:maxcolsize:N1;
        if(max(blocks1)~=N1), blocks1=[blocks1,N1]; end;
        blockcount1=numel(blocks1);
        
        blocks2=0:maxcolsize:N2;
        if(max(blocks2)~=N2), blocks2=[blocks2,N2]; end;
        blockcount2=numel(blocks2);
        
        D=zeros(N1,N2);
        
        for i=1:blockcount1-1
            blockrange1=(blocks1(i)+1):blocks1(i+1);
            for j=1:blockcount2-1
                blockrange2=(blocks2(j)+1):blocks2(j+1);
                
                D(blockrange1,blockrange2)=distDouble(varargin{1}(:,blockrange1),varargin{2}(:,blockrange2));
            end
        end
        
    end
		
end

function D = distSingle(singleMatrix)
K = singleMatrix'*singleMatrix;
vD = diag(K);
N = numel(vD);
A = repmat(vD(:),1,N);
D = -2*K + A + A';
end

function D = distDouble(firstMatrix,secondMatrix)
Dy = repmat(sum(secondMatrix.^2,1),size(firstMatrix,2),1);
Dx = repmat(sum(firstMatrix.^2,1)',1,size(secondMatrix,2));
Dxy = firstMatrix'*secondMatrix;
D = Dy + Dx - 2*Dxy;
end
