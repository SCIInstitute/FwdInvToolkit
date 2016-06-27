function InterpolatedCurve=InterpolateCurve(CurveParameters,InterpolationDensity,varargin)
% Author: Burak Erem

SizeCurveParameters=size(CurveParameters);
ProdSizeCurveParameters=prod(SizeCurveParameters);

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

% Form 1-D spline interpolation matrices of appropriate sizes
% TensorEdgeDimensions=sort(unique(SizeCurveParameters(2:end)));
% for i=TensorEdgeDimensions
%     SplineMatrix{i-2}=splineMatrix(i-2,1,InterpolationDensity).';
% end
TensorEdgeDimensions=SizeCurveParameters(2:end);
for i=1:numel(TensorEdgeDimensions)
    if(PeriodicityFlag(i)==0) % if not periodic
        SplineMatrix{i}=splineMatrix(TensorEdgeDimensions(i)-2,1,InterpolationDensity).';
    else % if periodic
        SplineMatrix{i}=splineMatrix(TensorEdgeDimensions(i)-2,1,InterpolationDensity,'Periodic').';
    end
end

% Initialize interpolated curves as curve parameters
InterpolatedCurve=CurveParameters;

% Interpolate spline curves using tensor-matrix "right" multiplication, for
% all the "right-hand" sides (i.e. tensor indices, not including the first one)
for i=1:numel(SizeCurveParameters)-1
    InterpolatedCurve=tensorRightMatrixMultiply(InterpolatedCurve,i,SplineMatrix{i});
end

