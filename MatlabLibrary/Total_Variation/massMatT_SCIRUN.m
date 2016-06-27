%% CALCULATES THE MASS MATRIX OF THE GIVEN GEOMETRY

% load data
	vertX = field1.node(:,1);
	vertY = field1.node(:,2);
	vertZ = field1.node(:,3);

% compute mass matrix
[ LocMassMat, jcbInt ] = BuildLocMassMat_LinearTet( rhsType,vertX, vertY, vertZ, ws, phi );
