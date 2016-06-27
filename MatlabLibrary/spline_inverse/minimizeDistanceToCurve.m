function [CurveParameters, Cost] = minimizeDistanceToCurve(InitialCurveParameters, DataPoints, InterpolationDensity, varargin)
% Expects InitialCurveParameters to be
% [DerivativeOfFirstKnot,KnotPointsInColumns,DerivativeOfLastKnot]
% TODO: more documentation
% Two optional arguments at end:
% If parameters end with ...,[mincostreduction],'JustDerivatives')
%    then the knot points will remain fixed, and the derivatives will be
%    changed during optimization.
% If the optional mincostreduction variable is specified, this will replace
%    the default minimum cost reduction (1e-6) for numerical convergence.
% 
% Author: Burak Erem
% Algorithm: The optimization algorithm used to solve the problem is
% similar in technique to the k-means algorithm (a.k.a. Lloyd's algorithm),
% with similar convergence properties (...yet to be formally proven).
% 

DEFAULTmincostreduction=1e-6;
justderivsflag=0;

SizeCurveParameters=size(InitialCurveParameters);
ProdSizeCurveParameters=prod(SizeCurveParameters);

PeriodicDims=[];
PeriodicityFlag=zeros(numel(SizeCurveParameters)-1,1);

% DataWeights=[];

if(numel(varargin)>0)
    if(~isempty(varargin{1}))
        mincostreduction=varargin{1};
    else
        mincostreduction=DEFAULTmincostreduction;
    end
    
    if(numel(varargin)>1)
        for i=2:numel(varargin)
            if(strcmpi(varargin{i},'JustDerivatives'))
                justderivsflag=1;
            end
            if(strcmpi(varargin{i},'Periodic'))
                % Unless this is followed by a numeric array of dimension
                % indices that should be periodic, assume they are all
                % periodic
                if(numel(varargin)>i)
                    if(~ischar(varargin{i+1}))
                        PeriodicDims=varargin{i+1};
                        PeriodicityFlag(PeriodicDims)=1;
                        if(numel(PeriodicityFlag)>numel(SizeCurveParameters)-1)
                            fprintf('WARNING: Dimensions specified as being periodic exceed input dimensions.\n')
                        end
                    end
                else
                    PeriodicityFlag=ones(size(PeriodicityFlag));
                end
            end
%             if(strcmpi(varargin{i},'Weights'))
%                 % Unless this is followed by a numeric array of weights,
%                 % give a warning and assume uniform weighting (the default)
%                 if(numel(varargin)>i)
%                     if(~ischar(varargin{i+1}))
%                         DataWeights=varargin{i+1};
%                         DataWeights=DataWeights(:).'; % make it a row vec
%                         if(numel(DataWeights)~=size(DataPoints,2))
%                             % TODO: make this a "real" Matlab error message
%                             fprintf('ERROR: Number of input data weights does not match number of input data points.\n')
%                             return
%                         end
%                     end
%                 else
%                     fprintf('WARNING: No data weights were provided! Assuming uniform weights (default setting)\n')
%                 end
%             end
        end
    end
else
    mincostreduction=DEFAULTmincostreduction;
end



% Minimize over CurveParameters -> S

% wolfeparams=struct([]);
wolfeparams.alphamax=1;

% phi=@(S)(TotalCurveProjections(reshape(S,SizeCurveParameters),DataPoints,InterpolationDensity,PeriodicDims,DataWeights));
phi=@(S)(TotalCurveProjections(reshape(S,SizeCurveParameters),DataPoints,InterpolationDensity,PeriodicDims));


% TensorEdgeDimensions=sort(unique(SizeCurveParameters(2:end)));
TensorEdgeDimensions=SizeCurveParameters(2:end);
for i=1:numel(TensorEdgeDimensions)
    if(PeriodicityFlag(i)==0) % if not periodic
        SplineMatrix{i}=splineMatrix(TensorEdgeDimensions(i)-2,1,InterpolationDensity).';
    else % if periodic
        SplineMatrix{i}=splineMatrix(TensorEdgeDimensions(i)-2,1,InterpolationDensity,'Periodic').';
    end
end

% dphidx=@(S)(full(reshape(CurveProjectionsDescentDirection(reshape(S,SizeCurveParameters),DataPoints,InterpolationDensity,SplineMatrix,justderivsflag,PeriodicDims,PeriodicityFlag,DataWeights),ProdSizeCurveParameters,1)));
dphidx=@(S)(full(reshape(CurveProjectionsDescentDirection(reshape(S,SizeCurveParameters),DataPoints,InterpolationDensity,SplineMatrix,justderivsflag,PeriodicDims,PeriodicityFlag),ProdSizeCurveParameters,1)));

[CurveParameters,dfdx,alpha]=steepestdescent(phi,dphidx,InitialCurveParameters(:),mincostreduction,wolfeparams);

CurveParameters=reshape(CurveParameters(:,end),SizeCurveParameters);
Cost=phi(CurveParameters);

end

% % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % %

% function Cost = TotalCurveProjections(CurveParameters,DataPoints,InterpolationDensity,Periodicity,DataWeights)
function Cost = TotalCurveProjections(CurveParameters,DataPoints,InterpolationDensity,Periodicity)
MAXELEMSPAIRWISEDISTS=5e6;

% 1. Interpolate the spline
if(sum(Periodicity(:))>0)
    ISet = InterpolateCurve(CurveParameters,InterpolationDensity,'Periodic',Periodicity);
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

% if(~isempty(DataWeights))
%     P2P=P2P.*repmat(DataWeights,size(P2P,1),1);
% end

% 4. Calculate total curve projection distances squared, divided by number of projected data points
Cost = (1/numel(DataPoints))*sum(min(P2P));

end

% % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % 

% function TargetChange = CurveProjectionsDescentDirection(CurveParameters,DataPoints,InterpolationDensity,SplineMatrix,justderivsflag,PeriodicDims,PeriodicityFlag,DataWeights)
function TargetChange = CurveProjectionsDescentDirection(CurveParameters,DataPoints,InterpolationDensity,SplineMatrix,justderivsflag,PeriodicDims,PeriodicityFlag)
SizeCurveParams=size(CurveParameters);
Ts=SizeCurveParams(2:end)-2;
TotalCurvePoints=InterpolationDensity.*(Ts-1)+Ts;
M=SizeCurveParams(1);

% Project DataPoints to the curves, returning the ProjectedPoints and the
% (linear) indices on the curves to which they project, ProjectedIndices
if(numel(PeriodicDims)>0)
    [ProjectedPoints,ProjectedIndices]=ProjectDataPointsToCurve(CurveParameters,DataPoints,InterpolationDensity,'Periodic',PeriodicDims);
else
    [ProjectedPoints,ProjectedIndices]=ProjectDataPointsToCurve(CurveParameters,DataPoints,InterpolationDensity);
end

% Compute a vector of desired changes, pointing from DataPoints to ProjectedPoints
DataPointsChange=ProjectedPoints-DataPoints;

% For every point on the curves, calculate the mean of the desired changes
TargetChange=zeros(M,prod(TotalCurvePoints));
for i=1:prod(TotalCurvePoints)
    if(sum(ProjectedIndices==i)>0)
        TargetChange(:,i)=mean(DataPointsChange(:,ProjectedIndices==i),2);
    end
end
% if(isempty(DataWeights)) % uniform average:
%     for i=1:prod(TotalCurvePoints)
%         if(sum(ProjectedIndices==i)>0)
%             TargetChange(:,i)=mean(DataPointsChange(:,ProjectedIndices==i),2);
%         end
%     end
% else % weighted average:
%     TargetWeights=zeros(prod(TotalCurvePoints),1);
%     for i=1:prod(TotalCurvePoints)
%         if(sum(ProjectedIndices==i)>0)
%             CurrWeights=repmat(DataWeights(:,ProjectedIndices==i),M,1)/sum(DataWeights(:,ProjectedIndices==i),2);
%             TargetChange(:,i)=sum(CurrWeights.*DataPointsChange(:,ProjectedIndices==i),2);
% %             TargetWeights(i)=sum(CurrWeights(1,:),2);
%             TargetWeights(i)=sum(sum(DataWeights(:,ProjectedIndices==i),2));
%         end
%     end
% end
    
% Reshape into a tensor
TargetChange=reshape(TargetChange,[M,TotalCurvePoints]);

% if(~isempty(DataWeights))
%     TargetWeights=reshape(TargetWeights,[1,TotalCurvePoints]);
%     TargetWeights=repmat(TargetWeights,[M,ones(size(TotalCurvePoints))]);
%     TargetWeights=TargetChange.*TargetWeights;
% end

% Find the least squares change in spline parameters for the mean desired changes
for i=1:numel(Ts)
    % is it periodic? if so, modify this so that we're not solving a
    % rank-deficient system
    if(PeriodicityFlag(i)==1) % Periodic?
        TargetChange=tensorRightMatrixSolve(TargetChange,i,SplineMatrix{i},1:size(SplineMatrix{i},1)-2);
    else % NOT Periodic!
        TargetChange=tensorRightMatrixSolve(TargetChange,i,SplineMatrix{i});
    end
    
%     if(isempty(DataWeights)) % uniform
%         % is it periodic? if so, modify this so that we're not solving a
%         % rank-deficient system
%         if(PeriodicityFlag(i)==1) % Periodic?
%             TargetChange=tensorRightMatrixSolve(TargetChange,i,SplineMatrix{i},1:size(SplineMatrix{i},1)-2);
%         else % NOT Periodic!
%             TargetChange=tensorRightMatrixSolve(TargetChange,i,SplineMatrix{i});
%         end
%     else % weighted
%         % is it periodic? if so, modify this so that we're not solving a
%         % rank-deficient system
%         RightHandSide=SplineMatrix{i}*diag(TargetWeights);
%         if(PeriodicityFlag(i)==1) % Periodic?
%             TargetChange=tensorRightMatrixSolve(TargetChange,i,RightHandSide,1:size(SplineMatrix{i},1)-2);
%         else % NOT Periodic!
%             TargetChange=tensorRightMatrixSolve(TargetChange,i,RightHandSide);
%         end
%     end

end

if(justderivsflag==1) % if this flag is enabled, we're going to generate a string with a command to be evaluated by matlab that zeros all knot point changes, leaving only derivative changes non-zero
    % NOTE: I'm not sure that zero-ing the search directions that we want
    % to remain unchanged is equivalent to solving for a least squares
    % TargetChange such that they are zero... will have to look into it.
    indexstring='TargetChange(:';
    for i=1:numel(Ts)
        indexstring=sprintf('%s,2:end-1',indexstring);
    end
    indexstring=sprintf('%s)=0;',indexstring);
    eval(indexstring);
end

end
