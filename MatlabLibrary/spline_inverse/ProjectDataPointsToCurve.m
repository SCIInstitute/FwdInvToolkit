function [ProjectedPoints,ProjectedIndices]=ProjectDataPointsToCurve(CurveParameters,DataPoints,InterpolationDensity,varargin)
% Author: Burak Erem

MAXELEMSPAIRWISEDISTS=5e6;

SizeCurveParameters=size(CurveParameters);
PeriodicityFlag=zeros(numel(SizeCurveParameters)-1,1);

if(numel(varargin)>0)
    for i=1:numel(varargin)
        if(strcmpi(varargin{i},'Periodic'))
            % Unless this is followed by a numeric array of dimension
            % indices that should be periodic, assume they are all
            % periodic
            if(numel(varargin)>i)
                if(~ischar(varargin{i+1}))
                    PeriodicityFlag(varargin{i+1})=1;
                    if(numel(PeriodicityFlag)>numel(SizeCurveParameters)-1)
                        fprintf('WARNING: Dimensions specified as being periodic exceed input dimensions.\n')
                    end
                end
            else
                PeriodicityFlag=ones(size(PeriodicityFlag));
            end
        end
    end
end

% 1. Interpolate the spline
if(sum(PeriodicityFlag(:))>0)
    ISet = InterpolateCurve(CurveParameters,InterpolationDensity,'Periodic',find(PeriodicityFlag));
else
    ISet = InterpolateCurve(CurveParameters,InterpolationDensity);
end

% 2. Reshape the tensor into a matrix with points in the columns
SizeISet = size(ISet);
ISet = reshape(ISet,SizeISet(1),prod(SizeISet(2:end)));

% 3. Calculate distances squared between ISet (row indices) and data (column indices)
if(size(ISet,2)*size(DataPoints,2)>MAXELEMSPAIRWISEDISTS)
    P2P = PairwiseDistance(ISet,DataPoints,floor(sqrt(MAXELEMSPAIRWISEDISTS)));
else
    P2P = PairwiseDistance(ISet,DataPoints);
end

% Find the index of the minimum of each column
[~,ProjectedIndices]=min(P2P);

% The projected points are the points that correspond to these minima
ProjectedPoints=ISet(:,ProjectedIndices);
