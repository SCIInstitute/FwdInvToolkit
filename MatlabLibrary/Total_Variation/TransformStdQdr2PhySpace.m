function [elmtX, elmtY, elmtZ] = TransformStdQdr2PhySpace(z_x, z_y, z_z, vertX, vertY, vertZ, option)
% this function maps the quadrature points from [-1,1] cubic domain to physical elements

nx = length( z_x);      ny = length(z_y);  nz = length( z_z);
elmtX = zeros( nx, ny, nz );  
elmtY = zeros( nx, ny, nz );  
elmtZ = zeros( nx, ny, nz );  

switch option
    case 'Tet'
        x1 = vertX(1);  x2 = vertX(2); x3 = vertX(3); x4 = vertX(4);   
        y1 = vertY(1);  y2 = vertY(2); y3 = vertY(3); y4 = vertY(4);   
        z1 = vertZ(1);   z2 = vertZ(2); z3 = vertZ(3);  z4 = vertZ(4);   
        for i = 1:nx
            cx = z_x( i );
            for j = 1: ny
                cy = z_y( j );
                for k = 1: nz
                    cz = z_z( k );
                    elmtX( i, j, k) = x1*(1-cx)*0.5 *(1-cy)*0.5*(1-cz)*0.5 + x2*(1+cx)*0.5*(1-cy)*0.5*(1-cz)*0.5...
                        + x3*(1+cy)*0.5*(1-cz)*0.5 + x4 * (1+cz)*0.5;
                    elmtY( i, j, k) = y1*(1-cx)*0.5 *(1-cy)*0.5*(1-cz)*0.5 + y2*(1+cx)*0.5*(1-cy)*0.5*(1-cz)*0.5...
                        + y3*(1+cy)*0.5*(1-cz)*0.5 + y4 * (1+cz)*0.5;
                    elmtZ( i, j, k) = z1*(1-cx)*0.5 *(1-cy)*0.5*(1-cz)*0.5 + z2*(1+cx)*0.5*(1-cy)*0.5*(1-cz)*0.5...
                        + z3*(1+cy)*0.5*(1-cz)*0.5 + z4 * (1+cz)*0.5;                    
                end
            end
        end
    case 'Prism' %the base isx-y plane, the height is z-axis 
        x1 = vertX(1);  x2 = vertX(2); x3 = vertX(3); x4 = vertX(4);    x5 = vertX(5);      x6 = vertX(6);
        y1 = vertY(1);  y2 = vertY(2); y3 = vertY(3); y4 = vertY(4);    y5 = vertY(5);      y6 = vertY(6);
        z1 = vertZ(1);   z2 = vertZ(2); z3 = vertZ(3);  z4 = vertZ(4);    z5 = vertZ(5);       z6 = vertZ(6);        
        for i = 1:nx
            cx = z_x(i);
            for j = 1:ny
                cy = z_y( j );
                for k = 1: nz
                    cz = z_z( k );
                    elmtX( i, j, k) = x1*(1-cx)*0.5 *(1-cy)*0.5*(1-cz)*0.5 + x2*(1+cx)*0.5*(1-cy)*0.5*(1-cz)*0.5...
                        + x3*(1+cy)*0.5*(1-cz)*0.5 + x4 *(1-cx)*0.5*(1-cy)*0.5*(1+cz)*0.5 + ...
                        x5*(1+cx)*0.5*(1-cy)*0.5*(1+cz)*0.5 + x6*(1+cy)*0.5*(1+cz)*0.5;            
                    
                    elmtY( i, j, k) = y1*(1-cx)*0.5 *(1-cy)*0.5*(1-cz)*0.5 + y2*(1+cx)*0.5*(1-cy)*0.5*(1-cz)*0.5...
                        + y3*(1+cy)*0.5*(1-cz)*0.5 + y4 *(1-cx)*0.5*(1-cy)*0.5*(1+cz)*0.5 + ...
                        y5*(1+cx)*0.5*(1-cy)*0.5*(1+cz)*0.5 + y6*(1+cy)*0.5*(1+cz)*0.5;              
                    
                    elmtZ( i, j, k) = z1*(1-cx)*0.5 *(1-cy)*0.5*(1-cz)*0.5 + z2*(1+cx)*0.5*(1-cy)*0.5*(1-cz)*0.5...
                        + z3*(1+cy)*0.5*(1-cz)*0.5 + z4 *(1-cx)*0.5*(1-cy)*0.5*(1+cz)*0.5 + ...
                        z5*(1+cx)*0.5*(1-cy)*0.5*(1+cz)*0.5 + z6*(1+cy)*0.5*(1+cz)*0.5;                         
                end
            end
        end
        
end
return;