function S = splineMatrix(NumberOfKnotPoints,DimensionOfSpace,InterpolationDensity,varargin)
% function S = splineMatrix(NumberOfKnotPoints,DimensionOfSpace,InterpolationDensity)
% Author: Burak Erem
% 
% Given a matrix X of [DerivAtFirstKnot,KnotPoints,DerivAtLastKnot],
% this function produces a matrix S that multiplies a vectorized X:
% y = S*vec(X)
% Y=reshape(y,DimensionOfSpace,T) is a matrix of points on the interpolated
% spline curve where
% T=InterpolationDensity*(NumberOfKnotPoints-1)+NumberOfKnotPoints
% is the number of interpolated points (i.e. with InterpolationDensity
% number of points between each knot point)
% 

PERIODIC=0;
if(numel(varargin)>0)
    if(strcmpi(varargin{1},'Periodic'))
        PERIODIC=1;
    end
end

N=NumberOfKnotPoints+2;

S=sparse(InterpolationDensity*(NumberOfKnotPoints-1)+NumberOfKnotPoints,N);

for i=1:N
    temp=zeros(1,N);
    temp(i)=1;
    
    ss=spline(1:NumberOfKnotPoints,temp);
    
    tempcol=ppval(ss,linspace(1,NumberOfKnotPoints,InterpolationDensity*(NumberOfKnotPoints-1)+NumberOfKnotPoints));
    
    S(:,i)=tempcol(:);
end

if(PERIODIC==1)
    % Make a matrix that will duplicate repeated values
    P=[eye(N-2,N);[0,1,zeros(1,N-2)];[1,zeros(1,N-1)]];
    S=S*P;
end

S=kron(S,speye(DimensionOfSpace));
