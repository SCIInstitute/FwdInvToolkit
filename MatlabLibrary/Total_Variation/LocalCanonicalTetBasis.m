function coefmat = LocalCanonicalTetBasis(option, order )
%Generate local expansion basis in canonical tetrahedron coord
%This function is largely outdated and replace by CalPhyBasisTet() and EvalBasisTet()
%Copyright: Dafang Wang, SCI Institute, 2008-01-09

%coefmat: each column is the coefficient for each basis polynomial
% option: 
%Base': order1: a1+ a2*x + a3*y + a4*z
%             order2: a1+ a2*x + a3*y + a4*z + a5*x^2 + a6*y^2 + a7*z^2 + a8* xy + a9*xz + a10*yz

%Deriv: the same order as 'Base
%


switch option
    case 'Base'
        switch order
            case 1
                NumLx = 4;
                Onode = [ 0,0,0; 1,0,0; 0,1,0; 0,0,1];
                Lb = [ ones(4,1), Onode];
                rhs = eye(NumLx);
                coefmat = Lb \ rhs;
            case 2
                NumLx = 10;
                Onode = [ 0 0 0;
                                 1 0 0;
                                 0 1 0;
                                 0 0 1;
                                 0.5 0 0;
                                 0 0.5 0;
                                 0 0 0.5;
                                 0.5 0.5 0;
                                 0.5 0 0.5;
                                 0 0.5 0.5;];
                Lb = zeros(NumLx, NumLx);             
                for i = 1: NumLx
                    x = Onode(i,1);  y=Onode(i,2);  z=Onode(i,3);
                    Lb(i,:) = [1, x, y, z, x*x, y*y, z*z, x*y, x*z, y*z];
                end
                rhs = eye(NumLx);
                coefmat = Lb \ rhs;
                
            case 3
        end
        
    case 'Derivative'
        switch order 
            case 1
                coefmat = cell{4,1};
                NumLx  = 4;
                coefmat{1:1:4} = zeros(4,3); %3 derivatives to x,y,z, 4 coefficients for linear basis
                coefmat{1} = [ [-1,0,0,0]', [-1,0,0,0]', [-1,0,0,0]'];
                coefmat{2} = [ [1,0,0,0]', [0,0,0,0]', [0,0,0,0]' ];
                coefmat{3} = [ [0,0,0,0]', [1,0,0,0]', [0,0,0,0]' ];
                coefmat{4} = [ [0,0,0,0]', [0,0,0,0]', [1,0,0,0]' ];
                
            case 2
            case 3
        end
end
    
return;