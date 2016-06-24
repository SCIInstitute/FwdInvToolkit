function coefmat = CalPhyBasisTet( order, vts )
% Each column of coefmat is the coefficients of the basis function defined
% in a prism in  physical domain. The vertices of the prism is given by vts 

% Input:  each row of vts gives the x/y/z coordinates of a vertice, the
% ordering of vertices must follow the right-hand rule:
% cross((v2-v1),(v3-v1)) = v4-v1

    switch order 
        case 1
            sz = 4;
            Lb = zeros(sz, sz);
                for i =1:sz
                    x = vts(i,1);
                    y = vts(i,2);
                    z = vts(i,3);
                    evec = [1, x, y, z];
                    Lb( i, :) = evec;
                end
             rhs = eye(sz);
             coefmat = Lb \ rhs;
             
        case 2
                sz = 10;
                Lb = zeros(sz, sz);
                Onode = zeros( sz, 3); %all vertices
                Onode(1:4,:) = vts;
                Onode( 5, :) = mean( vts([1,2], :) );
                Onode( 6, :) = mean( vts([1,3], :) );
                Onode( 7, :) = mean( vts([1,4], :) );
                Onode( 8, :) = mean( vts([2,3], :) );
                Onode( 9, :) = mean( vts([2,4], :) );
                Onode( 10, :) = mean( vts([3,4], :) );
                 for s = 1:sz
                    x = Onode(s,1);
                    y = Onode(s,2);
                    z = Onode(s,3);
                    evec = [ 1, x, y, z, x*x, y*y, z*z, x*y, x*z, y*z];
                    Lb( s, :) = evec;
                 end
                rhs = eye(sz);
                coefmat = Lb \ rhs;   
                
        case 3
%             disp('CalPhyBasisTet(): cubic case needed to be completed');
                sz = 20;
                Lb = zeros(sz, sz);
                Onode = zeros( sz, 3); %all vertices
                %node:
                Onode(1:4,:) = vts;
                %edge:
                Onode( 5, :) = (2*vts(1,:) + vts(2,:) ) / 3; %x^2
                Onode( 6, :) = (vts(1,:) + 2*vts(2,:) ) / 3; %x^3
                Onode( 7, :) =  (2*vts(1,:) + vts(3,:) ) / 3; %y^2
                Onode( 8, :) =  (vts(1,:) + 2*vts(3,:) ) / 3; %y^3
                Onode( 9, :) =  (2*vts(1,:) + vts(4,:) ) / 3; %z^2
                Onode( 10, :) =  (vts(1,:) + 2*vts(4,:) ) / 3; %z^3;
                
                Onode( 11, :) = (2*vts(2,:) + vts(3,:) ) /  3; %x^2 * y
                Onode( 12, :) = (vts(2,:) + 2*vts(3,:) ) /3; %x * y^2
                Onode( 13, :) =  (2*vts(2,:) + vts(4,:) )/3; %x^2 * z
                Onode( 14, :) =  (vts(2,:) + 2*vts(4,:) )/3; %x* z^2
                Onode( 15, :) =  (2*vts(3,:) + vts(4,:) ) /3; %y^2 * z
                Onode( 16, :) =  (vts(3,:) + 2*vts(4,:) ) /3; %y* z^2;
                
                %face:
                Onode(17, :) = mean( [Onode(1,:); Onode(2,:); Onode(3,:)] ); %xy
                Onode(18, :) = mean( [Onode(1,:); Onode(2,:); Onode(4,:)] );  %xz
                Onode(19, :) = mean( [Onode(1,:); Onode(3,:); Onode(4,:)] );  %yz
                Onode(20, :) = mean( [Onode(2,:); Onode(3,:); Onode(4,:)] ); %xYz
                
                 for s = 1:sz
                    x = Onode(s,1);
                    y = Onode(s,2);
                    z = Onode(s,3);
                    xx = x*x;   yy = y*y;  zz = z*z;
                    xy = x*y; yz=y*z; xz = x*z;
                    ev = [ 1, x, y, z, xx, xx*x, yy, yy*y,zz, zz*z, xx*y, x*yy, xx*z, x*zz, yy*z, y*zz, xy, xz, yz, x*y*z ];
                    Lb( s, :) = ev;
                 end
                rhs = eye(sz);
                coefmat = Lb \ rhs;   
    end
return