function [ LocMassMat, jcbInt ] = BuildLocMassMat_LinearTet( rhsType,vertX, vertY, vertZ, ws, phi )
% This function build the local mass matrix system. Only works for linear tet element.
% It hard codes the calculation of derivatives, quadratures are used only for integration
%Copyright: Dafang Wang, SCI Institute, 2008-01-30

%  rhsType: rright hand side type, 0 for laplace
%  elmtOrder  the order of the base test polynomials,1 for linear, 2 for quadratic
%  vertX: x-coordinates for 4 vertices of the tet
%  vertY, vertZ: same as vertX
%  zp, ws, the quadrature points and weights in x,y,z dimensions respectively. 
% 
%  IMPORTANT: 
% Basis func phi are values defined over quaduatures. It doesn't matter whether '-1to1' or '0to1' space. 
% The integration is done based on the '-1to1' tet. The weights have been set in that way.
%
%Output:

 w_x = ws(:,1);
  w_y = ws(:,2);
w_z = ws(:,3);
localMatDim = 4;
 
    %calculate the jacobi from [-1,1] std tet to physical tet elements
    [ jcbInt, jcbmatInt ] = JcbStdTet2GeneralTet( vertX, vertY, vertZ, '-1to1' );
    
    LocMassMat = zeros( localMatDim, localMatDim);
    
    for p = 1:localMatDim
        for q = p:localMatDim
            integrand = phi(:,:,:,p) .* phi(:,:,:,q);
             LocMassMat(p,q) = IntegrationInTet(integrand, w_x, w_y, w_z, jcbInt );
             LocMassMat(q,p) = LocMassMat(p,q);
        end
    end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the integration in the tet, follow the section 4.1.1.2 in spencer's book
% jcb, the jacobian from [-1,1] tet to physical tet, wx,wy,wz are weights
% derived from [-1,1] and have been adapted for tet integration
function r = IntegrationInTet(f, wx, wy, wz, jcb )
    nx = length(wx);  ny = length(wy);  nz = length(wz);
    for i = 1:nz
        f(:,:, i ) = f(:,:,i) * wz(i);
    end
    for i = 1:ny
        f(:,i,:) = f(:,i,:) * wy(i);
    end
    for i = 1:1:nx
        f(i,:,:) = f(i,:,:) * wx( i );
    end
    r = sum(sum(sum(f ) ) ) * abs(jcb );
return;