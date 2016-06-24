function feGradMag = ComputeGradMagPerTet( elmtH, feDerivMat, v)
%compute gradient magnitude for each linear tetrahedral element
% Input:
%   - v: the field data located at nodes
%   - feDerivMat:  a Nele*3*4 matrix, each 3*4 matric is the partial deriv for the 4 local basis functions in each ele
%       in 2D triangles, each tri has an 2*3 matrix for the partial deriv
eleNumH = size( elmtH, 1);
feGradMag = zeros( eleNumH, 1);
for i = 1:eleNumH
    cv = v(  elmtH(i,:) );
    feGradMag( i ) = norm( squeeze( feDerivMat(i,:,:) ) * cv, 2);
end

end