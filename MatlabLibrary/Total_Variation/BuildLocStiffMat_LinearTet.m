function [ LocStMat, LocDerivMat, vol, jcbInt, localRHS ] = BuildLocStiffMat_LinearTet( rhsFlag, ...
    conduct, vertX, vertY, vertZ, varargin)
% This function build the local stiffness matrix system
%   1. It hard codes the calculation of derivatives
%   2. The integration in the stiffness matrix is hard coded, as the derivatives are all zero.
%   3. If right-hand size source is needed, quadratures are used for its integration.
%   4. It differs from BuildLocStMatTet() in that it doesn't use quadrature for integration
%Copyright: Dafang Wang, SCI Institute, 2011-10-30

%  rhsType: right hand side type, 0 for laplace
%  rhsFlag:  1 if compute rhs
%  conduct:  the conductivity within this element
%  vertX: x-coordinates for 4 vertices of the tet
%  vertY, vertZ: same as vertX
% coefmatBase: each column is the coefficient of the basis in [0,1] canonical basis
% zp, ws, the quadrature points and weights in x,y,z dimensions respectively. 
%    Note that quadratures points are defined within a [-1,1] cube, but the
%    weights have been adapted for integration in tetrahedra

% Output:
%   - vol: element volume

localMatDim = 4;
LocStMat = zeros(4,4);        LocDerivMat = zeros(3,4);

%Calculate the coefficients of basis functions in the current tet, in physical space
coefmatBase = CalPhyBasisTet( 1, [vertX, vertY, vertZ] );      
for i = 1:4
    LocDerivMat(:, i ) = coefmatBase( 2:4, i ); %each column contains the dx, dy, dz for one basis func
end
%grad(u) = localDerivMat * [u1;u2;u3; u4]; the field value at 4 vertices
    
%calculate the jacobi from [-1,1] std tet to physical tet elements
[ jcbInt, jcbmatInt ] = JcbStdTet2GeneralTet( vertX, vertY, vertZ, '-1to1' );
vol = ComputeTetVolume( [vertX, vertY, vertZ] );

LocStMat = zeros( localMatDim, localMatDim);
for i =1:4
    for j = 1: i
        k = (conduct * LocDerivMat(:, j) )' * LocDerivMat(:, i)  * vol;
        LocStMat( i, j ) = k;
        LocStMat( j, i ) = k;
    end
end

if rhsFlag == 1 %need to compute rhs
    if ~exist('varargin')    disp('BuildLocStiffMat_LinearTet.m:  error \n');  return; end
    rhsType = varargin{1};
    if rhsType == 0
        localRHS = zeros( localMatDim, 1);
    else %use quadratuves
        zp = varargin{2};        ws = varargin{3};
        z_x = zp(:,1);      w_x = ws(:,1);
        z_y = zp(:,2);      w_y = ws(:,2);
        z_z = zp(:,3);      w_z = ws(:,3);
        qdrNumX = length(z_x);       qdrNumY = length(z_y);  qdrNumZ = length(z_z);
        elmtX = zeros(qdrNumX, qdrNumY, qdrNumZ);   
        [elmtX, elmtY, elmtZ] = TransformStdQdr2PhySpace(z_x, z_y, z_z, vertX, vertY, vertZ, 'Tet');
        vecXYZ = cat(4, elmtX, elmtY, elmtZ); %the last dimension is the xyz value

        [phi, phidX, phidY, phidZ] = EvalBasisTet(elmtOrder, coefmatBase, vecXYZ );
        
        %calculate the rhs source term
        eF = Test_Right_Side_Function( 3, rhsType, conduct, elmtX, elmtY, elmtZ );
        for p = 1:localMatDim
            rhsTerm = phi(:,:,:,p).* eF;
            localRHS( p,1 ) = IntegrationInTet(rhsTerm, w_x, w_y, w_z, jcbInt );
        end
    end
end

    
return;
end

function [vol] = ComputeTetVolume( vertices )

	vol = elemVolumeCalculation([1 2 3 4], vertices);
	
end
