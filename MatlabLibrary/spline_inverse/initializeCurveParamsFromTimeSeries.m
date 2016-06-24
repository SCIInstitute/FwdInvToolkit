function CurveParams = initializeCurveParamsFromTimeSeries(Data,NumberOfKnots)
% Author: Burak Erem

[N,T]=size(Data);

knotinds=floor(linspace(1,T,NumberOfKnots));
Knots=Data(:,knotinds);

FirstKnotDeriv=(Knots(:,2)-Knots(:,1))/(knotinds(2)-knotinds(1));
LastKnotDeriv=(Knots(:,end)-Knots(:,end-1))/(knotinds(end)-knotinds(end-1));

CurveParams=[FirstKnotDeriv,Knots,LastKnotDeriv];
