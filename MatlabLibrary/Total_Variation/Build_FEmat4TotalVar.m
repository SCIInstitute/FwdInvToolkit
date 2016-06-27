%Build FE-based matrices to be used for total variation
clear;
filename = '../ts_ht1_mesh1_Tr1.mat';
FileOut = 'tvFE_mesh1_Tr1.mat';
load(filename);

feDerivMat = zeros( eleNumH, 3, 4); % partial deriv in each element
feStiffMat = zeros( eleNumH, 4, 4); %stiffness matrix in each cell

% stH = sparse( ndNumH, ndNumH);
conduct = eye(3,3);
eleVol_H = zeros( eleNumH,1);

%%%%%%%%%%%%% feDerivMat and feStiffMat

%-----------------Build the stiffness matrix. -----------------
for i = 1:eleNumH
    ch = elmtH( i, 1:4);
    vx = ndMatH( ch, 1);   vy = ndMatH(ch, 2);    vz = ndMatH(ch, 3);
  
    [ LocStMat, LocDerivMat, vol ] = BuildLocStiffMat_LinearTet( 0, conduct, vx, vy, vz);
    eleVol_H( i ) = vol;
    feDerivMat( i, :, :) = LocDerivMat;
    feStiffMat( i, :, :) = LocStMat;
%     stH( ch, ch) = stH(ch, ch) + LocStMat;
end
% full( MatrixDiff( stH, stMatH_std)) %validation

%-----------------Build the mass matrix. -----------------
[z_x, w_x] = JacobiGLZW( 4, 0, 0);
[z_y, w_y] = JacobiGRZW( 4, 1, 0 );     w_y = w_y / 2;
[z_z, w_z] =  JacobiGRZW( 4, 2, 0 );     w_z = w_z / 4;
zpTet = [ z_x, z_y, z_z];          wsTet = [w_x, w_y, w_z];
qdTet = Transform2StdTetSpace( z_x, z_y, z_z, '0to1');
coefmatBaseTet = LocalCanonicalTetBasis('Base', 1);
phiTet = EvalBasisTet(1, coefmatBaseTet, qdTet);

gMassMat = sparse( ndNum, ndNum);
for i = 1:eleNum
    cn = Elmt( i, :);
    vx = NdMat(cn, 1);  vy = NdMat(cn,2);       vz = NdMat(cn, 3);
     [ locMassMat, jcbInt ] = BuildLocMassMat_LinearTet( 0,vx, vy, vz, wsTet, phiTet );
     gMassMat( cn, cn) = gMassMat(cn, cn) + locMassMat;
end

%-----------------Build the FE matrix for dual -----------------
[z_x, w_x] = JacobiGLZW( 4, 0, 0);
[z_y, w_y] = JacobiGRZW( 4, 1, 0 );     w_y = w_y / 2;
[z_z, w_z] =  JacobiGRZW( 4, 2, 0 );     w_z = w_z / 4;
zpTet = [ z_x, z_y, z_z];          wsTet = [w_x, w_y, w_z];
qdTet = Transform2StdTetSpace( z_x, z_y, z_z, '0to1');
coefmatBaseTet = LocalCanonicalTetBasis('Base', 1);
phiTet = EvalBasisTet(1, coefmatBaseTet, qdTet);

Ax = sparse( ndNumH, ndNumH);       Ay = Ax;        Az = Ax;
for i = 1:eleNumH
    ch = elmtH( i, 1:4);
    vx = ndMatH( ch, 1);   vy = ndMatH(ch, 2);    vz = ndMatH(ch, 3);
    locDerivMat = squeeze(feDerivMat( i, :, :));
    [ locAx, locAy, locAz ] = BuildLocCrossDerivMat_LinearTet( vx, vy, vz, wsTet, phiTet, locDerivMat );
    
    Ax( ch, ch) = Ax( ch, ch) + locAx;
    Ay( ch, ch) = Ay( ch, ch) + locAy;
    Az( ch, ch) = Az( ch, ch) + locAz;
end

massMatT = full(gMassMat( ndv_TsSurf, ndv_TsSurf));
M1 = mR' * (gStMat \ mQ'); %M1 is full matrix
M1A = M1 \ [ Ax, Ay, Az];     

%-----------------------  End -------------------------------------------------
clear i z_* w_* cn ch Loc*Mat loc*Mat vx vy vz opt jcb* locA* vol
save( FileOut );

% %Compute true TV value
% feGradMag = zeros( eleNumH, 1);  %gradient magnitude
% feGradMag0 = ComputeGradMagPerTet( elmtH, feDerivMat, vTMP0); %true grad mag
% feTV0 = dot( eleVol_H, feGradMag0);



