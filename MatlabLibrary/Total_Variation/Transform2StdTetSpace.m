function VecXYZ = Transform2StdTetSpace( z_x, z_y, z_z, option)
% Transform quadrature points from [-1,1] * [-1,1] to std tet space
% z_x, z_y, z_z : 1D quadrature points in x,y,z coordinates
% VecXYZ(:,:,1) for x coord, VecXYZ(:,:,2) for y coord, VecXYZ(:,:,3) for z coord

nx = length( z_x);      ny = length(z_y);  nz = length( z_z);
VecXYZ = zeros( nx, ny, nz, 3);

switch option
    case '0to1'
            for i = 1: nx
                cx = z_x( i );
                for j = 1: ny
                    cy = z_y( j );
                    for k = 1:nz
                        cz = z_z( k );
                        VecXYZ( i, j, k, 1) = (1+cx)*0.5 * (1-cy) * 0.5 * ( 1-cz) * 0.5;
                        VecXYZ( i, j, k, 2) = (1+cy) * 0.5 * ( 1-cz) * 0.5;
                        VecXYZ( i, j, k, 3) = (1+cz) * 0.5;
                    end
                end
            end
            
    case '-1to1'
end
return;