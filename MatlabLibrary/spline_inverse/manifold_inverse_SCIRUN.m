%% HELP: 
%
%		This code implements the inverse solutions pipeline presented in
%		the paper:
%			Erem, Coll-font, Martinez Orellana - 2013 - 
%			Using Transmural Regularization and Dynamic Modeling 
%			for Non-Invasive Cardiac Potential Imaging of 
%			Endocardial Pacing Sites with Imprecise Thoracic Geometry.
%
%		INPUTS:
%			- i1 - <N,T>double - measured potentials on the torso.
%			- i2 - <N,M>double - forward matrix.
%			- i3 - <L,M>double - regularization matrix.
%			- i4 - <3,1>int - regularization constant params.
%								i4(1) - log10 min lambda.
%								i4(2) - log10 max lambda
%								i4(3) - num lambda.
%
%		OUTPUT:
%			- o1 - <M,T>double - estimated heart potentials.
%
%		AUTHOR:
%			Jaume Coll-Font <jcollfont@gmail.com>
%			based on the original code by Burak Erem
%
%


%% LOAD
	
	torsoPotentials = i1;							% load torso potentials
	A = i2;											% load forward matrx
	R = i3;											% load regularization matrix
	regparamrange=10.^linspace(i4(1),i4(2),1e4);	% select lambda range: 10^i4(1)...10^i4(2)

%% PARAMETERS
	
	% signal pre-filtering
		lowpasswin=10;							% low pass window

	% curve interpolation
		NumberOfKnots=12;						% number of knots
		InterpolationDensity=100;				% interpolation density
		ProjectionInterpolationDensity=200;		% final projection density
		minderivcostreduction=500;				% minimum reconstruct err for deriv
		minoverallcostreduction= 1e-6;			% minimum reconstr err total

%% PRE_FILTERING
	torsoPotentials=lowpassma(torsoPotentials,lowpasswin);
	
%% FIT 1D MANIFOLD ON TORSO DATA
	% Solve for initial derivatives
	fprintf('Fitting first and last curve derivatives...\n')
	CurveParams = initializeCurveParamsFromTimeSeries(torsoPotentials,NumberOfKnots);
	CurveParams = minimizeDistanceToCurve(CurveParams,torsoPotentials,InterpolationDensity,minderivcostreduction,'JustDerivatives');

	% Solve for all parameters (find knot points)
	fprintf('Fitting all curve parameters...\n')
	CurveParams = minimizeDistanceToCurve(CurveParams,torsoPotentials,InterpolationDensity,minoverallcostreduction);

	% Project the data to the curve
	fprintf('Applying Tikhonov inverse to curve parameters...\n')
	[ torsoPotentials_Manifold,  torsoPotentials_timeWarp ] = ProjectDataPointsToCurve(CurveParams,torsoPotentials,ProjectionInterpolationDensity);
		
%% APPLY INVERSE TIKHONOV TO KNOT POINTS
	CurveParams_heart = tikhonovburak(A,R,CurveParams,regparamrange);
	
%% RECONSTRUCT HEART POTENTIALS FROM MANIFOLD
	heartPotentials_Manifold = InterpolateCurve( CurveParams_heart, ProjectionInterpolationDensity);
	heartPotentials_reconstructed = heartPotentials_Manifold( :, torsoPotentials_timeWarp); 
		
%% OUTPUT
	o1 = heartPotentials_reconstructed;