function [phi, phidX, phidY, phidZ] = EvalBasisTet(baseOrder, coefmat, VecXYZ)
%this function evaluates the derivatives of basis functions in an arbitrary tetrahedral domain
%Input:
% -coefmat: each column is the coefficients of one basis function, the organization of coefficient depends on 
%  the polynomial order of basis function
% -VecXYZ is a 4D array containing the quadrature points which have been
%   transformed into the current tet domain
%Copyright: Dafang Wang, SCI Institute, 2008-11-21

    nx = size(VecXYZ,1);
    ny = size(VecXYZ,2);
    nz = size(VecXYZ,3);
    
    basisNum = size( coefmat, 2);
    phidX = zeros( nx, ny, nz, basisNum);     phidY=phidX;        phidZ = phidX;      phi = phidX;
    
    if nargout == 1 %just evaluate function, no derivative
        for s = 1:basisNum
            coef = coefmat(:, s);
            for i = 1:nx
                for j = 1:ny
                    for k = 1:nz
                        cx = VecXYZ(i,j,k,1);       cy = VecXYZ(i,j,k,2);       cz = VecXYZ(i,j,k,3);
                        phi(i,j,k,s) = subEvalBasisTet(baseOrder, coef, [cx;cy;cz] );
                    end
                end
            end
        end
    else
        for s = 1:basisNum
            coef = coefmat(:, s);
            for i = 1:nx
                for j = 1:ny
                    for k = 1:nz
                        cx = VecXYZ(i,j,k,1);       cy = VecXYZ(i,j,k,2);       cz = VecXYZ(i,j,k,3);
                        phi(i,j,k,s) = subEvalBasisTet(baseOrder, coef, [cx;cy;cz] );
                        dvec =  subEvalBasisDerivTet(baseOrder, coef, [cx;cy;cz] );
                        phidX(i,j,k,s) = dvec(1);
                        phidY(i,j,k,s) = dvec(2);
                        phidZ(i,j,k,s) = dvec(3);
                    end
                end
            end
        end
    end
    
return;

function vec = subEvalBasisTet(bsOdr, bcth, coord)
%this function evaluate the basis function, bcth is the coefficients
    x = coord(1);
    y = coord(2);
    z = coord(3);
    
    switch bsOdr
      case 1
        e = [1; x; y; z];
      case 2
        e = [1; x; y; z; x*x; y*y; z*z; x*y; x*z; y*z];
      case 3
         xx = x*x;  yy=y*y;  zz = z*z;
         xy = x*y;  xz = x*z;  yz = y*z;
         e = [1; x; y; z; xx; xx*x; yy; yy*y; zz; zz*z; xx*y; x*yy; xx*z; x*zz; yy*z; y*zz; xy; xz; yz; x*y*z]; 
    end

    vec = dot(bcth, e);
return


function dvec = subEvalBasisDerivTet(bsOdr, bcth, coord)
    x = coord(1);
    y = coord(2);
    z = coord(3);
    switch bsOdr
      case 1
        e1 = [0; 1; 0; 0];
        e2 = [0; 0; 1; 0];
        e3 = [0; 0; 0; 1];
      case 2
        e1 = [0; 1; 0; 0; 2*x; 0; 0; y; z; 0]; %dx
        e2 = [0; 0; 1; 0; 0; 2*y; 0; x; 0; z]; %dy
        e3 = [0; 0; 0; 1; 0; 0; 2*z; 0; x; y]; %dz
      case 3
        e1 = [0; 1; 0; 0; 2*x; 3*x*x;  0; 0; 0; 0; 2*y*x; y*y; 2*z*x; z*z; 0; 0; y; z; 0; y*z ]; %dx
        e2 = [0; 0; 1; 0; 0; 0; 2*y; 3*y*y; 0; 0; x*x; x*2*y; 0; 0; 2*z*y; z*z; x; 0; z; x*z]; %dy
        e3 = [0; 0; 0; 1; 0; 0; 0; 0; 2*z; 3*z*z; 0; 0; x*x; 2*x*z; y*y; 2*y*z; 0; x; y; x*y]; %dz
     
    end
    dvec = zeros(3,1);
    dvec(1) = dot(bcth, e1);
    dvec(2) = dot(bcth, e2);
    dvec(3) = dot(bcth, e3);
return